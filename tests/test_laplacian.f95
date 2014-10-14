include 'testutil.i'

program test_laplacian
  use adjac
  use testutil
  implicit none

  integer, parameter :: n = 5

  double precision, dimension(n) :: y_value
  double precision, dimension(n,n) :: jac

  type(adjac_double), dimension(n) :: x
  type(adjac_double), dimension(n) :: y
  integer i, j

  do j = 1, n
     call adjac_set_independent(x(j), dble(j), j)
  end do

  ! Some calculation, for example the Laplacian plus nonlinearity
  y(1) = x(1)
  y(n) = x(n)
  do j = 2, n-1
     y(j) = x(j-1) - 2d0*x(j) + x(j+1)
  end do

  call adjac_get_value(y, y_value)
  write(*,*) y_value

  call adjac_get_dense_jacobian(y, jac)
  do i = 1, n
     write(*,*) jac(i,:)
  end do
end program test_laplacian
