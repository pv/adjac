program main
  use adjac

  integer, parameter :: n = 100, filliter = 10
  
  double precision :: xval(n), yval(n), J(n,n)
  integer :: rep, i

  do i = 1, n
     xval(i) = 1d0 + i
  end do

  do rep = 1, int(1e8/(n*n))
     call doit(xval, yval, J)
  end do

  do i = 1, 5
     write(*,*) yval(i)
  end do
  do i = 1, 5
     write(*,*) real(J(i,1:5))
  end do

contains
  subroutine doit(xval, yval, J)
    implicit none
    double precision, intent(in) :: xval(n)
    double precision, intent(out) :: yval(n), J(n,n)
    double precision :: x(n), y(n)
    integer :: i
    double precision, parameter :: eps = 1d-7

    call laplacian(xval, yval)

    x = xval

    do i = 1, n
        x(i) = xval(i) + eps
        call laplacian(x, y)
        J(:,i) = (y - yval) / eps
        x(i) = xval(i)
    end do
  end subroutine doit

  subroutine laplacian(x, y)
    implicit none
    double precision, dimension(n), intent(in) :: x
    double precision, dimension(n), intent(out) :: y
    integer :: i, p

    y(1) = -2*x(1) + x(2)
    y(n) = -2*x(n) + x(n-1)
    do i = 2, n-1
       y(i) = x(i+1) - 2*x(i) + x(i-1)
    end do

    do p = 1, filliter
       do i = 2, n-1
          y(i) = y(i-1) - 2*y(i) + y(i+1)/(3.1415d0 + y(i))
       end do
    end do

  end subroutine laplacian
end program main
