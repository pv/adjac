include 'testutil.i'

program test_laplacian
  use adjac
  use testutil
  implicit none

  integer, parameter :: n = 5

  double precision, dimension(n) :: y_value
  double precision, dimension(n,n) :: jac
  double precision, allocatable, dimension(:) :: jac_val
  integer, allocatable, dimension(:) :: jac_i, jac_j

  double precision, dimension(n) :: p, xval, dy, dy2
  type(adjac_double), dimension(n) :: x
  type(adjac_double), dimension(n) :: y
  integer i, j

  call adjac_reset()
  do j = 1, n
     call adjac_set_independent(x(j), dble(j), j)
  end do

  ! Some calculation, for example the Laplacian plus nonlinearity
  call laplacian(x, y)
  call adjac_get_value(y, y_value)
  write(*,*) y_value

  ! Dense jacobian
  call adjac_get_dense_jacobian(y, jac)
  do i = 1, n
     write(*,*) jac(i,:)
  end do

  ! Obtain jacobian in sparse coordinate (i, j, value) format
  write(*,*) '-- COO laplacian'
  call adjac_get_coo_jacobian(y, jac_val, jac_i, jac_j)
  do j = 1, size(jac_i,1)
     write(*,*) jac_i(j), jac_j(j), real(jac_val(j))
  end do

  ! Evaluate jacobian-vector products
  do j = 1, n
     p(j) = 1d0/(1 + j)
     xval(j) = dble(j)
  end do
  call adjac_reset(.true.)
  call adjac_set_independent(x, xval, p)
  call laplacian(x, y)
  call adjac_get_value(y, y_value, dy)
  dy2 = matmul(jac, p)

  if (maxval(abs(dy - dy2)) > 1d-12) then
     write(*,*) 'jac product FAIL'
  else
     write(*,*) 'jac product OK'
  end if
contains
  subroutine laplacian(x, y)
    implicit none 
    type(adjac_double), dimension(:), intent(in) :: x
    type(adjac_double), dimension(:), intent(out) :: y
    integer :: n
    n = size(x)
    y(1) = x(1)
    y(n) = x(n)
    do j = 2, n-1
       y(j) = x(j-1) - 2d0*x(j) + x(j+1)
    end do
  end subroutine laplacian
end program test_laplacian
