! =====
! adjac
! =====
!
! Automatic Differentiation (forward-mode) for sparse Jacobians.
!
! Usage::
!
!    subroutine somefunction(x_value, y_value, jac_dense, jac_i, jac_j, jac_val)
!        use adjac
!        implicit none
!    
!        double precision, dimension(n), intent(in) :: x_value
!        double precision, dimension(n), intent(out) :: y_value
!        double precision, dimension(n, n), intent(out) :: jac_dense
!        integer, allocatable, dimension(:), intent(inout) :: jac_i, jac_j, jac_ind, jac_indptr
!        double precision, allocatable, dimension(:), intent(inout) :: jac_val, jac_val2
!    
!        type(adjac_double), dimension(n) :: x
!        type(adjac_double), dimension(n) :: y
!        integer :: j, nnz
!    
!        do j = 1, n
!            call adjac_set_independent(x(j), x_value(j), j)
!        end do
!    
!        ! Some calculation, for example the Laplacian plus nonlinearity
!        y(1) = x(1)
!        y(n) = x(n)
!        do j = 2, n-1
!           y(j) = x(j-1) - 2*x(j) + x(j+1) + cos(x(j))
!        end do
!
!        ! Obtain function value
!        call adjac_get_value(y, y_value)
!
!        ! Obtain jacobian in dense format
!        call adjac_get_dense_jacobian(y, jac_dense)
!
!        ! Obtain jacobian in sparse coordinate (i, j, value) format
!        nnz = adjac_get_nnz(y)
!        allocate(jac_i(nnz), jac_j(nnz), jac_val(nnz))
!        call adjac_get_coo_jacobian(y, jac_i, jac_j, jac_val)
!
!        ! Obtain jacobian in compressed sparse row format
!        allocate(jac_ind(n+1), jac_indptr(nnz), jac_val2(nnz))
!        call adjac_get_csr_jacobian(y, jac_ind, jac_indptr, jac_val2)
!    end subroutine somefunction

! Copyright (c) 2014, Pauli Virtanen <pav@iki.fi>
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its
! contributors may be used to endorse or promote products derived from
! this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

module adjac
  private

  type, public :: adjac_double
     double precision :: value
     integer :: n
     integer, allocatable, dimension(:) :: taylor_i
     double precision, allocatable, dimension(:) :: taylor_val
  end type adjac_double

  type, public :: adjac_complex_an
     double complex :: value
     integer :: n
     integer, allocatable, dimension(:) :: taylor_i
     double complex, allocatable, dimension(:) :: taylor_val
  end type adjac_complex_an

  type, public :: adjac_complex
     type(adjac_double) :: re, im
  end type adjac_complex

  interface adjac_get_value
     module procedure get_value_one, get_value_many
  end interface adjac_get_value

  public assignment(=)
  interface assignment(=)
     module procedure assign_ad, assign_bd, assign_bz, assign_ba
  end interface

  public operator(+)
  interface operator(+)
     module procedure add_aa, add_ad, add_da
     module procedure add_bb, add_bz, add_zb
     module procedure add_ab, add_ba
     module procedure add_az, add_za
     module procedure add_bd, add_db
  end interface

  public operator(-)
  interface operator(-)
     module procedure sub_aa, sub_ad, sub_da
     module procedure sub_bb, sub_bz, sub_zb
     module procedure sub_ab, sub_ba
     module procedure sub_az, sub_za
     module procedure sub_bd, sub_db
     module procedure neg_a, neg_b
  end interface

  public operator(*)
  interface operator(*)
     module procedure mul_aa, mul_ad, mul_da
     module procedure mul_bb, mul_bz, mul_zb
     module procedure mul_ab, mul_ba
     module procedure mul_az, mul_za
     module procedure mul_bd, mul_db
  end interface operator(*)

  public operator(/)
  interface operator(/)
     module procedure div_aa, div_ad, div_da
     module procedure div_bb, div_bz, div_zb
     module procedure div_ab, div_ba
     module procedure div_az, div_za
     module procedure div_bd, div_db
  end interface operator(/)

  public matmul
  interface matmul
     module procedure matmul_aa, matmul_ad, matmul_da
     module procedure matmul_bb, matmul_bz, matmul_zb
  end interface matmul

  public dble
  interface dble
     module procedure dble_a, dble_b
  end interface dble

  public aimag
  interface aimag
     module procedure aimag_b
  end interface aimag

  public conjg
  interface conjg
     module procedure conjg_b
  end interface conjg

  public exp
  interface exp
     module procedure exp_a, exp_b
  end interface exp

  public sin
  interface sin
     module procedure sin_a, sin_b
  end interface sin

  public cos
  interface cos
     module procedure cos_a, cos_b
  end interface cos

  public log
  interface log
     module procedure log_a, log_b
  end interface log

  public adjac_set_independent, adjac_get_value, &
       adjac_get_dense_jacobian, adjac_get_csr_jacobian, &
       adjac_get_coo_jacobian, adjac_get_nnz

contains
  subroutine adjac_set_independent(x, xval, j)
    implicit none
    type(adjac_double), intent(inout) :: x
    double precision, intent(in) :: xval
    integer, intent(in) :: j

    allocate(x%taylor_i(1))
    allocate(x%taylor_val(1))

    x%n = 1
    x%value = xval
    x%taylor_i(1) = j
    x%taylor_val(1) = 1
  end subroutine adjac_set_independent

  subroutine get_value_one(y, val)
    implicit none
    type(adjac_double), intent(in) :: y
    double precision, intent(out) :: val
    val = y%value
  end subroutine get_value_one

  subroutine get_value_many(y, val)
    implicit none
    type(adjac_double), dimension(:), intent(in) :: y
    double precision, dimension(size(y,1)), intent(out) :: val
    integer :: j
    do j = 1, size(val,1)
       val(j) = y(j)%value
    end do
  end subroutine get_value_many

  function adjac_get_nnz(y) result(nnz)
    type(adjac_double), dimension(:), intent(in) :: y
    integer :: nnz
    nnz = 0
    do i = 1, size(y,1)
       nnz = nnz + y(i)%n
    end do
  end function adjac_get_nnz

  subroutine adjac_get_dense_jacobian(y, jac_dense)
    implicit none
    type(adjac_double), dimension(:), intent(inout) :: y
    double precision, dimension(:,:), intent(out) :: jac_dense

    integer :: i, p, j

    jac_dense = 0

    do i = 1, size(y,1)
       do p = 1, y(i)%n
          jac_dense(i, y(i)%taylor_i(p)) = jac_dense(i, y(i)%taylor_i(p)) + y(i)%taylor_val(p)
       end do
    end do
  end subroutine adjac_get_dense_jacobian

  subroutine adjac_get_coo_jacobian(y, jac_i, jac_j, jac_val)
    implicit none
    type(adjac_double), dimension(:), intent(inout) :: y
    integer, dimension(:), intent(out) :: jac_i, jac_j
    double precision, dimension(:), intent(out) :: jac_val

    integer :: i, p, j, k

    k = 1
    do i = 1, size(y,1)
       if (y(i)%n > 0) then
          jac_i(k:k+y(i)%n-1) = i
          jac_j(k:k+y(i)%n-1) = y(i)%taylor_i(1:y(i)%n)
          jac_val(k:k+y(i)%n-1) = y(i)%taylor_val(1:y(i)%n)
          k = k + y(i)%n
       end if
    end do
  end subroutine adjac_get_coo_jacobian

  subroutine adjac_get_csr_jacobian(y, jac_ind, jac_indptr, jac_val)
    implicit none
    type(adjac_double), dimension(:), intent(inout) :: y
    integer, dimension(:), intent(out) :: jac_ind, jac_indptr
    double precision, dimension(:), intent(out) :: jac_val

    integer :: i, p, j, k

    k = 1
    jac_ind(1) = 1
    do i = 1, size(y,1)
       if (y(i)%n > 0) then
          jac_indptr(k:k+y(i)%n-1) = y(i)%taylor_i(1:y(i)%n)
          jac_val(k:k+y(i)%n-1) = y(i)%taylor_val(1:y(i)%n)
          k = k + y(i)%n
       end if
       jac_ind(i+1) = k
    end do
  end subroutine adjac_get_csr_jacobian

  pure subroutine axpy_sparse_vectors(a, na, ia, va, nb, ib, vb, nc, ic, vc)
    ! c := a*x + y
    implicit none
    double precision, intent(in) :: a
    integer, intent(in) :: na, nb
    integer, intent(inout) :: nc
    integer, dimension(na), intent(in) :: ia
    integer, dimension(nb), intent(in) :: ib
    integer, dimension(nc), intent(out) :: ic
    double precision, dimension(na), intent(in) :: va
    double precision, dimension(nb), intent(in) :: vb
    double precision, dimension(nc), intent(out) :: vc

    double precision :: tmp
    integer :: p, q, r

    ! Sum sorted, until either va or vb exhausted
    p = 1
    q = 1
    r = 1
    do
       if (p > na .or. q > nb) exit

       if (ia(p) < ib(q)) then
          vc(r) = a * va(p)
          ic(r) = ia(p)
          p = p + 1
          r = r + 1
       else if (ia(p) > ib(q)) then
          vc(r) = vb(q)
          ic(r) = ib(q)
          q = q + 1
          r = r + 1
       else
          tmp = a * va(p) + vb(q)
          if (tmp .ne. 0) then
             vc(r) = tmp
             ic(r) = ia(p)
             r = r + 1
          end if
          p = p + 1
          q = q + 1
       end if
    end do

    ! Copy rest
    do p = p, na, 1
       vc(r) = a * va(p)
       ic(r) = ia(p)
       r = r + 1
    end do
    do p = q, nb, 1
       vc(r) = vb(q)
       ic(r) = ib(q)
       r = r + 1
    end do

    ! Finish
    nc = r - 1
  end subroutine axpy_sparse_vectors

  !--------------------------------------------------------------------------
  ! Overloaded operators
  !--------------------------------------------------------------------------

  !!
  !! assignment(=)
  !!

  pure elemental subroutine assign_ad(x, y)
    implicit none
    type(adjac_double), intent(out) :: x
    double precision, intent(in) :: y
    x%value = y
    x%n = 0
  end subroutine assign_ad

  pure elemental subroutine assign_bd(x, y)
    implicit none
    type(adjac_complex), intent(out) :: x
    double precision, intent(in) :: y
    x%re = y
    x%im = 0d0
  end subroutine assign_bd

  pure elemental subroutine assign_bz(x, y)
    implicit none
    type(adjac_complex), intent(out) :: x
    double complex, intent(in) :: y
    x%re = dble(y)
    x%im = aimag(y)
  end subroutine assign_bz

  pure elemental subroutine assign_ba(x, y)
    implicit none
    type(adjac_complex), intent(out) :: x
    type(adjac_double), intent(in) :: y
    x%re = y
    x%im = 0d0
  end subroutine assign_ba

  !!
  !! operator(+)
  !!

  ! X + Y = x + y + (x_j - y_j) dj

  pure elemental function add_aa(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x, y
    type(adjac_double) :: z

    z%value = x%value + y%value

    z%n = x%n + y%n
    allocate(z%taylor_i(z%n))
    allocate(z%taylor_val(z%n))

    call axpy_sparse_vectors(&
         1d0, &
         x%n, x%taylor_i, x%taylor_val, &
         y%n, y%taylor_i, y%taylor_val, &
         z%n, z%taylor_i, z%taylor_val)
  end function add_aa

  pure elemental function add_ad(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_double) :: z
    z = x
    z%value = z%value + y
  end function add_ad

  pure elemental function add_da(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_double) :: z
    z = y + x
  end function add_da

  pure elemental function add_bb(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re + y%re
    z%im = x%im + y%im
  end function add_bb

  pure elemental function add_bz(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re + dble(y)
    z%im = x%im + aimag(y)
  end function add_bz

  pure elemental function add_zb(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) + y%re
    z%im = aimag(x) + y%im
  end function add_zb

  pure elemental function add_ba(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re + y
    z%im = x%im
  end function add_ba

  pure elemental function add_ab(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x + y%re
    z%im = y%im
  end function add_ab

  pure elemental function add_az(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x + dble(y)
    z%im = aimag(y)
  end function add_az

  pure elemental function add_za(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) + y
    z%im = aimag(x)
  end function add_za

  pure elemental function add_bd(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_complex) :: z
    z = x + dcmplx(y)
  end function add_bd

  pure elemental function add_db(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = dcmplx(x) + y
  end function add_db

  !!
  !! operator(-)
  !!

  ! X - Y = x - y + (x_j - y_j) dj

  pure elemental function sub_aa(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x, y
    type(adjac_double) :: z

    integer :: n

    z%value = x%value - y%value

    z%n = x%n + y%n
    allocate(z%taylor_i(z%n))
    allocate(z%taylor_val(z%n))

    call axpy_sparse_vectors(&
         -1d0, &
         y%n, y%taylor_i, y%taylor_val, &
         x%n, x%taylor_i, x%taylor_val, &
         z%n, z%taylor_i, z%taylor_val)
  end function sub_aa

  pure elemental function sub_ad(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_double) :: z
    z = x
    z%value = z%value - y
  end function sub_ad

  pure elemental function sub_da(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_double) :: z
    z = y
    z%value = x - z%value
    z%taylor_val = -z%taylor_val
  end function sub_da

  pure elemental function sub_bb(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re - y%re
    z%im = x%im - y%im
  end function sub_bb

  pure elemental function sub_bz(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re - dble(y)
    z%im = x%im - aimag(y)
  end function sub_bz

  pure elemental function sub_zb(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) - y%re
    z%im = aimag(x) - y%im
  end function sub_zb

  pure elemental function sub_ba(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re - y
    z%im = x%im
  end function sub_ba

  pure elemental function sub_ab(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x - y%re
    z%im = -y%im
  end function sub_ab

  pure elemental function sub_az(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x - dble(y)
    z%im = -aimag(y)
  end function sub_az

  pure elemental function sub_za(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) - y
    z%im = aimag(x)
  end function sub_za

  pure elemental function sub_bd(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_complex) :: z
    z = x - dcmplx(y)
  end function sub_bd

  pure elemental function sub_db(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = dcmplx(x) - y
  end function sub_db

  !!
  !! operator(-), unary
  !!

  pure elemental function neg_a(x) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_double) :: z
    z = 0d0 - x
  end function neg_a

  pure elemental function neg_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex) :: z
    z = (0d0,0d0) - z
  end function neg_b

  !!
  !! operator(*)
  !!

  ! X*Y = x*y + (x y_j + y x_j) dj

  pure elemental function mul_aa(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x, y
    type(adjac_double) :: z
    z = (x * y%value) + (x%value * y)
    z%value = x%value * y%value
  end function mul_aa

  pure elemental function mul_ad(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_double) :: z
    if (y == 0) then
       z%n = 0
       z%value = 0
    else
       z = x
       z%value = z%value * y
       z%taylor_val(1:z%n) = z%taylor_val(1:z%n) * y
    end if
  end function mul_ad

  pure elemental function mul_da(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_double) :: z
    z = mul_ad(y, x)
  end function mul_da

  pure elemental function mul_bb(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re * y%re - x%im * y%im
    z%im = x%re * y%im + x%im * y%re
  end function mul_bb

  pure elemental function mul_bz(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re * dble(y) - x%im * aimag(y)
    z%im = x%re * aimag(y) + x%im * dble(y)
  end function mul_bz

  pure elemental function mul_zb(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) * y%re - aimag(x) * y%im
    z%im = dble(x) * y%im + aimag(x) * y%re
  end function mul_zb

  pure elemental function mul_ba(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re * y
    z%im = x%im * y
  end function mul_ba

  pure elemental function mul_ab(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x * y%re
    z%im = x * y%im
  end function mul_ab

  pure elemental function mul_az(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x * dble(y)
    z%im = x * aimag(y)
  end function mul_az

  pure elemental function mul_za(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) * y
    z%im = aimag(x) * y
  end function mul_za

  pure elemental function mul_bd(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_complex) :: z
    z = x * dcmplx(y)
  end function mul_bd

  pure elemental function mul_db(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = dcmplx(x) * y
  end function mul_db

  !!
  !! operator(/)
  !!

  ! X/Y = x/y + (x_j/y - x y_j/y**2) dj

  pure elemental function div_aa(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x, y
    type(adjac_double) :: z
    z = ((1d0 / y%value) * x) + ((-x%value / (y%value**2)) * y)
    z%value = x%value / y%value
  end function div_aa

  pure elemental function div_ad(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_double) :: z
    z = (1d0 / y) * x
  end function div_ad

  pure elemental function div_da(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_double) :: z
    z = (-x / (y%value**2)) * y
    z%value = x / y%value
  end function div_da

  pure elemental function div_bb(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = x * conjg(y) / (dble(y)*dble(y) + aimag(y)*aimag(y))
  end function div_bb

  pure elemental function div_bz(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z = x * conjg(y) / (dble(y)*dble(y) + aimag(y)*aimag(y))
  end function div_bz

  pure elemental function div_zb(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = x * conjg(y) / (dble(y)*dble(y) + aimag(y)*aimag(y))
  end function div_zb

  pure elemental function div_ba(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re / y
    z%im = x%im / y
  end function div_ba

  pure elemental function div_ab(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = x * conjg(y) / (dble(y)*dble(y) + aimag(y)*aimag(y))
  end function div_ab

  pure elemental function div_az(x, y) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    double complex, intent(in) :: y
    type(adjac_complex) :: z
    z = x * conjg(y) / (dble(y)*dble(y) + aimag(y)*aimag(y))
  end function div_az

  pure elemental function div_za(x, y) result(z)
    implicit none
    double complex, intent(in) :: x
    type(adjac_double), intent(in) :: y
    type(adjac_complex) :: z
    z%re = dble(x) / y
    z%im = aimag(x) / y
  end function div_za

  pure elemental function div_bd(x, y) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    double precision, intent(in) :: y
    type(adjac_complex) :: z
    z%re = x%re / y
    z%im = x%im / y
  end function div_bd

  pure elemental function div_db(x, y) result(z)
    implicit none
    double precision, intent(in) :: x
    type(adjac_complex), intent(in) :: y
    type(adjac_complex) :: z
    z = dcmplx(x) / y
  end function div_db

  !!
  !! matmul
  !!

  function matmul_aa(x, y) result(z)
    implicit none
    type(adjac_double), dimension(:,:), intent(in) :: x, y
    type(adjac_double), dimension(size(x,1),size(y,2)) :: z

    integer i, j, k

    if (size(x,2) .ne. size(y,1)) then
       write(*,*) 'invalid array sizes in matmul'
       stop
    end if

    do j = 1, size(y,2)
       do i = 1, size(x,1)
          z(i,j) = x(i,1)*y(1,j)
          do k = 2, size(x,2)
             z(i,j) = z(i,j) + x(i,k)*y(k,j)
          end do
       end do
    end do
  end function matmul_aa

  function matmul_ad(x, y) result(z)
    implicit none
    type(adjac_double), dimension(:,:), intent(in) :: x
    double precision, dimension(:,:), intent(in) :: y
    type(adjac_double), dimension(size(x,1),size(y,2)) :: z

    integer i, j, k

    if (size(x,2) .ne. size(y,1)) then
       write(*,*) 'invalid array sizes in matmul'
       stop
    end if

    do j = 1, size(y,2)
       do i = 1, size(x,1)
          z(i,j) = x(i,1)*y(1,j)
          do k = 2, size(x,2)
             z(i,j) = z(i,j) + x(i,k)*y(k,j)
          end do
       end do
    end do
  end function matmul_ad

  function matmul_da(x, y) result(z)
    implicit none
    double precision, dimension(:,:), intent(in) :: x
    type(adjac_double), dimension(:,:), intent(in) :: y
    type(adjac_double), dimension(size(x,1),size(y,2)) :: z

    integer i, j, k

    if (size(x,2) .ne. size(y,1)) then
       write(*,*) 'invalid array sizes in matmul'
       stop
    end if

    do j = 1, size(y,2)
       do i = 1, size(x,1)
          z(i,j) = x(i,1)*y(1,j)
          do k = 2, size(x,2)
             z(i,j) = z(i,j) + x(i,k)*y(k,j)
          end do
       end do
    end do
  end function matmul_da

  function matmul_bb(x, y) result(z)
    implicit none
    type(adjac_complex), dimension(:,:), intent(in) :: x, y
    type(adjac_complex), dimension(size(x,1),size(y,2)) :: z

    integer i, j, k

    if (size(x,2) .ne. size(y,1)) then
       write(*,*) 'invalid array sizes in matmul'
       stop
    end if

    do j = 1, size(y,2)
       do i = 1, size(x,1)
          z(i,j) = x(i,1)*y(1,j)
          do k = 2, size(x,2)
             z(i,j) = z(i,j) + x(i,k)*y(k,j)
          end do
       end do
    end do
  end function matmul_bb

  function matmul_bz(x, y) result(z)
    implicit none
    type(adjac_complex), dimension(:,:), intent(in) :: x
    double complex, dimension(:,:), intent(in) :: y
    type(adjac_complex), dimension(size(x,1),size(y,2)) :: z

    integer i, j, k

    if (size(x,2) .ne. size(y,1)) then
       write(*,*) 'invalid array sizes in matmul'
       stop
    end if

    do j = 1, size(y,2)
       do i = 1, size(x,1)
          z(i,j) = x(i,1)*y(1,j)
          do k = 2, size(x,2)
             z(i,j) = z(i,j) + x(i,k)*y(k,j)
          end do
       end do
    end do
  end function matmul_bz

  function matmul_zb(x, y) result(z)
    implicit none
    double complex, dimension(:,:), intent(in) :: x
    type(adjac_complex), dimension(:,:), intent(in) :: y
    type(adjac_complex), dimension(size(x,1),size(y,2)) :: z

    integer i, j, k

    if (size(x,2) .ne. size(y,1)) then
       write(*,*) 'invalid array sizes in matmul'
       stop
    end if

    do j = 1, size(y,2)
       do i = 1, size(x,1)
          z(i,j) = x(i,1)*y(1,j)
          do k = 2, size(x,2)
             z(i,j) = z(i,j) + x(i,k)*y(k,j)
          end do
       end do
    end do
  end function matmul_zb

  !!
  !! dble
  !!

  pure elemental function dble_a(x) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_double) :: z
    z = x
  end function dble_a

  pure elemental function dble_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_double) :: z
    z = x%re
  end function dble_b

  !!
  !! aimag
  !!

  pure elemental function aimag_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_double) :: z
    z = x%im
  end function aimag_b

  !!
  !! conjg
  !!

  pure elemental function conjg_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex) :: z
    z%re = x%re
    z%im = -x%im
  end function conjg_b

  !!
  !! exp
  !!

  pure elemental function exp_a(x) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_double) :: z
    double precision :: v, dv
    v = exp(z%value)
    dv = v
    z = dv*x
    z%value = v
  end function exp_a

  pure elemental function exp_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex) :: z
    double complex :: v, dv
    v = exp(dcmplx(x%re%value, x%im%value))
    dv = v
    z = dv*x
    z%re%value = dble(v)
    z%im%value = aimag(v)
  end function exp_b

  !!
  !! sin
  !!

  pure elemental function sin_a(x) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_double) :: z
    double precision :: v, dv
    v = sin(x%value)
    dv = cos(x%value)
    z = dv*x
    z%value = v
  end function sin_a

  pure elemental function sin_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex) :: z
    double complex :: v, dv
    v = sin(dcmplx(x%re%value, x%im%value))
    dv = cos(dcmplx(x%re%value, x%im%value))
    z = dv*x
    z%re%value = dble(v)
    z%im%value = aimag(v)
  end function sin_b

  !!
  !! cos
  !!

  pure elemental function cos_a(x) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_double) :: z
    double precision :: v, dv
    v = cos(x%value)
    dv = -sin(x%value)
    z = dv*x
    z%value = v
  end function cos_a

  pure elemental function cos_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex) :: z
    double complex :: v, dv
    v = cos(dcmplx(x%re%value, x%im%value))
    dv = -sin(dcmplx(x%re%value, x%im%value))
    z = dv*x
    z%re%value = dble(v)
    z%im%value = aimag(v)
  end function cos_b

  !!
  !! log
  !!

  pure elemental function log_a(x) result(z)
    implicit none
    type(adjac_double), intent(in) :: x
    type(adjac_double) :: z
    double precision :: v, dv
    v = log(x%value)
    dv = 1d0/x%value
    z = dv*x
    z%value = v
  end function log_a

  pure elemental function log_b(x) result(z)
    implicit none
    type(adjac_complex), intent(in) :: x
    type(adjac_complex) :: z
    double complex :: v, dv
    v = log(dcmplx(x%re%value, x%im%value))
    dv = 1d0/dcmplx(x%re%value, x%im%value)
    z = dv*x
    z%re%value = dble(v)
    z%im%value = aimag(v)
  end function log_b
end module adjac
