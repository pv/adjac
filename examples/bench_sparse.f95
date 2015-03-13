program main
  use adjac

  integer, parameter :: blocksize = 128, n = blocksize*1000, filliter = 3

  double precision :: xval(n), yval(n)
  double precision, allocatable, dimension(:) :: jac_val
  integer, allocatable, dimension(:) :: jac_indices, jac_indptr
  integer :: rep, i

  do i = 1, n
     xval(i) = 1d0 + i
  end do

  do rep = 1, 5
     if (rep > 1) then
        deallocate(jac_val, jac_indices, jac_indptr)
     end if
     call doit(xval, yval, jac_val, jac_indices, jac_indptr)
  end do

  write(*,*) 'nnz =', size(jac_val)
  write(*,*) 'sparsity =', (size(jac_val)*1d0)/(1d0*n)/(1d0*n)
  write(*,*) jac_val(:5)
  write(*,*) jac_indices(:5)
  write(*,*) jac_indptr(:5)

contains
  subroutine doit(xval, yval, jac_val, jac_indices, jac_indptr)
    implicit none
    double precision, intent(in) :: xval(n)
    double precision, intent(out) :: yval(n)
    type(adjac_double) :: x(n), y(n)
    integer :: i, nnz
    double precision, allocatable, dimension(:) :: jac_val
    integer, allocatable, dimension(:) :: jac_indices, jac_indptr

    call adjac_reset()
    call adjac_set_independent(x, xval)
    call oper(x, y)

    call adjac_get_value(y, yval)
    nnz = adjac_get_nnz(y)
    allocate(jac_val(nnz), jac_indices(nnz), jac_indptr(n+1))
    call adjac_get_csr_jacobian(y, jac_val, jac_indices, jac_indptr)
  end subroutine doit

  subroutine oper(x, y)
    implicit none
    type(adjac_double), dimension(n), intent(in) :: x
    type(adjac_double), dimension(n), intent(out) :: y
    integer :: i, k, p

    y(1) = -2*x(1) + x(2)
    y(n) = -2*x(n) + x(n-1)
    do i = 2, n-1
       y(i) = x(i+1) - 2*x(i) + x(i-1)
    end do

    do k = 1, n, blocksize
       do p = 1, filliter
          do i = k+1, k + blocksize-1
             y(i) = y(i-1) - 2*y(i) + y(i+1)/(3.1415d0 + y(i))
          end do
       end do
    end do
  end subroutine oper
end program main
