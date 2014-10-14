    subroutine my_func(x, y)
        implicit none
        double complex, dimension(3), intent(in) :: x
        double complex, dimension(2), intent(out) :: y

        integer :: j
        do j = 1, 2
            y(j) = log(x(j) / (1 + cos(x(j+1))**2))
        end do
    end subroutine my_func

    subroutine my_func_jac(x_value, y_value, dy_dx)
        use adjac
        implicit none
        double complex, dimension(3), intent(in) :: x_value
        double complex, dimension(2), intent(out) :: y_value
        double complex, dimension(2,3), intent(out) :: dy_dx

	type(adjac_complexan), dimension(3) :: x
	type(adjac_complexan), dimension(2) :: y
        integer :: j

	call adjac_set_independent(x, x_value)

        do j = 1, 2
            y(j) = log(x(j) / (1 + cos(x(j+1))**2))
        end do

	call adjac_get_value(y, y_value)
	call adjac_get_dense_jacobian(y, dy_dx)
    end subroutine my_func_jac

program simple
  implicit none
  double precision, parameter :: dx = 1e-7

  double complex, dimension(3) :: x, x2
  double complex, dimension(2) :: y, y2
  double complex, dimension(2,3) :: dy_dx

  integer :: i

  do i = 1, 3
     x(i) = i
  end do

  ! Evaluate function values
  call my_func(x, y)
  write(*,*) '-- my_func:'
  write(*,*) y

  ! Find jacobian by numerical differentiation
  dy_dx = 0
  do i = 1, 3
     x2 = x
     x2(i) = x(i) + dx
     call my_func(x2, y2)
     dy_dx(:,i) = (y2 - y) / dx
  end do

  write(*,*) '-- Jacobian via numerical differentiation:'
  do i = 1, 2
     write(*,*) cmplx(dy_dx(i,:))
  end do

  ! Get the same results via adjac
  dy_dx = 0
  call my_func_jac(x, y, dy_dx)
  write(*,*) '-- my_func_jac value:'
  write(*,*) y
  write(*,*) '-- my_func_jac Jacobian:'
  do i = 1, 2
     write(*,*) cmplx(dy_dx(i,:))
  end do

end program simple
