program main
  use adjac

  integer, parameter :: blocksize = 128, n = blocksize*1000, filliter = 3

  double precision :: xval(n), yval(n)
  double precision, allocatable, dimension(:) :: jac_val
  integer, allocatable, dimension(:) :: jac_i, jac_j
  integer :: rep, i

  do i = 1, n
     xval(i) = 1d0 + i
  end do

  do rep = 1, 5
     if (rep > 1) then
        deallocate(jac_val, jac_i, jac_j)
     end if
     call doit(xval, yval, jac_val, jac_i, jac_j)
  end do

  write(*,*) 'nnz =', size(jac_val)
  write(*,*) 'sparsity =', (size(jac_val)*1d0)/(1d0*n)/(1d0*n)
  write(*,*) jac_val(1:5)
  write(*,*) jac_i(1:5)
  write(*,*) jac_j(1:5)

contains
  subroutine doit(xval, yval, jac_val, jac_i, jac_j)
    implicit none
    double precision, intent(in) :: xval(n)
    double precision, intent(out) :: yval(n)
    type(adjac_double) :: x(n), y(n)
    double precision, allocatable, dimension(:), intent(inout) :: jac_val
    integer, allocatable, dimension(:), intent(inout) :: jac_i, jac_j
    integer :: nnz

    call adjac_reset()
    call adjac_set_independent(x, xval)
    call oper(x, y)

    call adjac_get_value(y, yval)
    call adjac_get_coo_jacobian(y, nnz, jac_val, jac_i, jac_j)
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
