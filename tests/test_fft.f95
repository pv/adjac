program test_fft
  use adjac
  use adjac_fft
  implicit none

  integer, parameter :: n = 128
  type(adjac_double), dimension(2*n) :: x
  double precision, dimension(2*n) :: x_value
  double precision, dimension(2*n,2*n) :: jac
  double complex :: v1, v2
  integer :: i, j

  do i = 1, n
     x_value(2*i-1) = 1 + i + i**2
     x_value(2*i) = 0
     call adjac_set_independent(x(2*i-1), x_value(2*i-1), 2*i-1)
     call adjac_set_independent(x(2*i), x_value(2*i), 2*i)
  end do

  write(*,*) 'forward transform'
  call fft(n, x_value)
  call fft(n, x)
  call adjac_get_dense_jacobian(x, jac)

  do i = 1, n
     v1 = x(i)%value
     v2 = x_value(i)
     if (abs(v1 - v2) > 1d-10) then
        write(*,*) 'FAIL', i, v1, v2
        stop
     end if
     do j = 1, n
        v1 = jac(2*i-1,2*j-1) + (0,1)*jac(2*i,2*j-1)
        v2 = exp(-2*(0,3.141592653589793d0)*(i-1)*(j-1)/n)
        if (abs(v1 - v2) > 1d-10) then
           write(*,*) 'FAIL', i, j, v1, v2
           stop
        end if
     end do
  end do
  write(*,*) 'OK'

  write(*,*) 'inverse transform'
  call ifft(n, x_value)
  call ifft(n, x)
  call adjac_get_dense_jacobian(x, jac)

  do i = 1, n
     if (mod(i,2) == 1) then
        v2 = 1 + (1 + i/2) + (1 + i/2)**2
     else
        v2 = 0
     end if

     v1 = x(i)%value
     if (abs(v1 - v2) > 1d-10) then
        write(*,*) 'FAIL', 'avalue', i, v1, v2
        stop
     end if

     v2 = x_value(i)
     if (abs(v1 - v2) > 1d-10) then
        write(*,*) 'FAIL', 'value', i, v1, v2
        stop
     end if

     do j = 1, n
        v1 = jac(2*i-1,2*j-1) + (0,1)*jac(2*i,2*j-1)
        if (i == j) then
           v2 = 1
        else
           v2 = 0
        end if
        if (abs(v1 - v2) > 1d-10) then
           write(*,*) 'FAIL', i, j, v1, v2
           stop
        end if
     end do
  end do
  write(*,*) 'OK'
end program test_fft
