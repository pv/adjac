include 'testutil.i'

program test_sqrt
  implicit none
  write(*,*) '-- a'
  call test_sqrt_a()
  write(*,*) '-- b'
  call test_sqrt_b()
end program test_sqrt

subroutine test_sqrt_a()
  use adjac
  use testutil
  implicit none

  double precision :: x_value, res_value(1), res_dvalue(1,1)
  type(adjac_double) :: x(1), res(1)

  call adjac_reset()

  x_value = -4d0
  call adjac_set_independent(x(1), x_value, 1)

  res = abs(x)

  call adjac_get_value(res, res_value)
  call adjac_get_dense_jacobian(res, res_dvalue)

  write(*,*) 'res =', res_value
  call assert_equal(res_value(1), abs(x_value))

  write(*,*) 'jac =', res_dvalue
  call assert_equal(res_dvalue(1,1), -1d0)
end subroutine test_sqrt_a

subroutine test_sqrt_b()
  use adjac
  use testutil
  implicit none

  double complex :: x_value
  double precision :: res_value(1), res_dvalue(1,2)
  type(adjac_double) :: x_parts(2)
  type(adjac_complex) :: x(1)
  type(adjac_double) :: res(1)

  call adjac_reset()

  x_value = (3d0, 4d0)
  
  call adjac_set_independent(x_parts(1), dble(x_value), 1)
  call adjac_set_independent(x_parts(2), aimag(x_value), 2)

  x(1)%re = x_parts(1)
  x(1)%im = x_parts(2)

  res = abs(x)

  call adjac_get_value(res, res_value)
  call adjac_get_dense_jacobian(res, res_dvalue)

  write(*,*) 'res =', res_value
  call assert_equal(res_value(1), abs(x_value))

  write(*,*) 'jac =', res_dvalue
  call assert_equal(res_dvalue(1,1), dble(x_value) / abs(x_value))
  call assert_equal(res_dvalue(1,2), aimag(x_value) / abs(x_value))
end subroutine test_sqrt_b
