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
  call assert_equal(y_value(1), dble(1))
  do j = 2, n-1
     call assert_equal(y_value(j), dble(j-1 - 2*j + j+1))
  end do
  call assert_equal(y_value(n), dble(n))

  call adjac_get_dense_jacobian(y, jac)
  do i = 1, n
     do j = 1, n
        write(*,*) i, j
        if (i == j) then
           if (i == 1 .or. i == n) then
              call assert_equal(jac(i,j), 1d0)
           else
              call assert_equal(jac(i,j), -2d0)
           end if
        else if (i .ne. 1 .and. i .ne. n .and. (i == j+1 .or. i == j-1)) then
           call assert_equal(jac(i,j), 1d0)
        else
           call assert_equal(jac(i,j), 0d0)
        end if
     end do
  end do
end program test_laplacian
