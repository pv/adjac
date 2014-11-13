program bench_advection
  implicit none
  integer, parameter :: nx = 100
  integer, parameter :: nt = 100
  integer, parameter :: nr = 500

  double precision, parameter :: pi = 4.0*atan(1.0)
  double precision :: q_init(nx)
  
  double precision, parameter :: dt = 0.125
  double precision :: jacobian(nx,nx)

  integer :: i

  do i = 1, nx
     q_init(i) = (0.5d0+0.5d0*sin(((i-1)*2d0*pi)/(nx-1.5d0)))+0.0001d0
  end do

  do i = 1, nr
     call doit(q_init, jacobian)
     if (i == 1) then
        write(*,*) jacobian(1,1:5)
        write(*,*) jacobian(2,1:5)
        write(*,*) jacobian(3,1:5)
        write(*,*) jacobian(4,1:5)
        write(*,*) jacobian(5,1:5)
     end if
  end do

contains

  subroutine doit(q_init_values, jacobian)
    use adjac
    implicit none
    double precision, intent(in) :: q_init_values(nx)
    double precision, intent(out) :: jacobian(nx,nx)

    type(adjac_double) :: q_init(nx), q(nx)

    call adjac_reset()
    call adjac_set_independent(q_init, q_init_values)

    call toon(nt, dt, q_init, q)
    call adjac_get_dense_jacobian(q, jacobian)
  end subroutine doit

  subroutine lax_wendroff(nt, c, q_init, q)
    use adjac
    implicit none
    integer, intent(in) :: nt
    double precision, intent(in) :: c
    type(adjac_double), intent(in) :: q_init(nx)
    type(adjac_double), intent(inout) :: q(nx)

    type(adjac_double) :: flux(nx-1)
    integer :: i, j

    do i = 1, nx
       q(i) = q_init(i)
    end do

    do j = 1, nt
       do i = 1, nx-1
          flux(i) = 0.5d0*c*(q(i)+q(i+1)+c*(q(i)-q(i+1)))
       end do
       do i = 2, nx-1
          q(i) = q(i) + flux(i-1) - flux(i)
       end do
       q(1) = q(nx-1)
       q(nx) = q(2)
    end do
  end subroutine lax_wendroff

  subroutine toon(nt, c, q_init, q)
    use adjac
    implicit none
    integer, intent(in) :: nt
    double precision, intent(in) :: c
    type(adjac_double), intent(in) :: q_init(nx)
    type(adjac_double), intent(inout) :: q(nx)

    type(adjac_double) :: flux(nx-1)
    integer :: i, j

    do i = 1, nx
       q(i) = q_init(i)
    end do

    do j = 1, nt
       do i = 1, nx-1
          flux(i) = (exp(c*log(q(i)/q(i+1)))-1.0d0) * q(i)*q(i+1) / (q(i)-q(i+1))
       end do
       do i = 2, nx-1
          q(i) = q(i) + flux(i-1) - flux(i)
       end do
       q(1) = q(nx-1)
       q(nx) = q(2)
    end do
  end subroutine toon
end program bench_advection
