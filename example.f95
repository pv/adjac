program example
  use adjac
  implicit none
  
  integer, parameter :: n = 3*2
 
  double precision, dimension(n) :: x_value
  double precision, dimension(n) :: y_value
  integer, allocatable, dimension(:) :: jac_i, jac_j
  double precision, allocatable, dimension(:) :: jac_val
  
  type(adjac_double), dimension(n) :: x
  type(adjac_double), dimension(n) :: y
  type(adjac_complex), dimension(3) :: yx
  integer j, nnz
  
  do j = 1, n
     x_value(j) = 1d0
     call adjac_set_independent(x(j), x_value(j), j)
  end do

  ! Some calculation, for example the Laplacian plus nonlinearity
  y(1) = x(1)
  y(n) = x(n)
  do j = 2, n-1
     y(j) = x(j-1) - 2d0*x(j) + x(j+1) + cos(x(j))
  end do

  ! Obtain function value
  write(*,*) '--'
  call adjac_get_value(y, y_value)
  write(*,*) real(y_value)
  
  ! Obtain jacobian in sparse coordinate (i, j, value) format
  write(*,*) '--'
  nnz = adjac_get_nnz(y)
  allocate(jac_i(nnz), jac_j(nnz), jac_val(nnz))
  call adjac_get_coo_jacobian(y, jac_i, jac_j, jac_val)
  do j = 1, size(jac_i,1)
     write(*,*) jac_i(j), jac_j(j), real(jac_val(j))
  end do
end program example
