include 'testutil.i'

program test_laplacian
  use adjac
  use testutil
  implicit none

  integer, parameter :: n = 5

  double precision, dimension(n) :: y_value
  double precision, dimension(n,n) :: jac
  double precision, allocatable, dimension(:) :: jac_val, jac_val2
  integer, allocatable, dimension(:) :: jac_i, jac_j, jac_indices, jac_indptr

  type(adjac_double), dimension(n) :: x
  type(adjac_double), dimension(n) :: y
  integer i, j, nnz

  call adjac_reset()
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

  ! Dense jacobian
  call adjac_get_dense_jacobian(y, jac)
  do i = 1, n
     write(*,*) jac(i,:)
  end do

  ! Obtain jacobian in sparse coordinate (i, j, value) format
  write(*,*) '-- COO laplacian'
  nnz = adjac_get_nnz(y)
  allocate(jac_val(nnz), jac_i(nnz), jac_j(nnz))
  call adjac_get_coo_jacobian(y, jac_val, jac_i, jac_j)
  do j = 1, size(jac_i,1)
     write(*,*) jac_i(j), jac_j(j), real(jac_val(j))
  end do

  ! Obtain jacobian in sparse CSR format
  write(*,*) '-- CSR laplacian'
  allocate(jac_val2(nnz), jac_indices(nnz), jac_indptr(n+1))
  call adjac_get_csr_jacobian(y, jac_val2, jac_indices, jac_indptr)
  write(*,*) jac_val2
  write(*,*) jac_indices
  write(*,*) jac_indptr
end program test_laplacian
