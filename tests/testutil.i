! -*-f90-*-
module testutil
  double precision, parameter :: default_eps = 1d-8

  interface assert_equal
     module procedure assert_equal_i, assert_equal_d, assert_equal_z
  end interface

  interface assert_close
     module procedure assert_close_i, assert_close_d, assert_close_z
  end interface

contains

  subroutine assert_equal_i(a, b)
    implicit none
    integer, intent(in) :: a, b
    if (a.ne.b) then
       write(*,*) 'FAIL', a, '!=', b
    end if
  end subroutine assert_equal_i

  subroutine assert_equal_d(a, b)
    implicit none
    double precision, intent(in) :: a, b
    if (a.ne.b) then
       write(*,*) 'FAIL', a, '!=', b
    end if
  end subroutine assert_equal_d

  subroutine assert_equal_z(a, b)
    implicit none
    complex(kind=kind(0d0)), intent(in) :: a, b
    if (a.ne.b) then
       write(*,*) 'FAIL', a, '!=', b
    end if
  end subroutine assert_equal_z

  subroutine assert_close_i(a, b)
    implicit none
    integer, intent(in) :: a, b
    double precision :: eps = default_eps
    if (abs(a-b).gt.eps) then
       write(*,*) 'FAIL', a, '!=', b, '+-', eps
    end if
  end subroutine assert_close_i

  subroutine assert_close_d(a, b)
    implicit none
    double precision, intent(in) :: a, b
    double precision :: eps = default_eps
    if (abs(a-b).gt.eps) then
       write(*,*) 'FAIL', a, '!=', b, '+-', eps
    end if
  end subroutine assert_close_d

  subroutine assert_close_z(a, b)
    implicit none
    complex(kind=kind(0d0)), intent(in) :: a, b
    double precision :: eps = default_eps
    if (abs(a-b).gt.eps) then
       write(*,*) 'FAIL', a, '!=', b, '+-', eps
    end if
  end subroutine assert_close_z
end module testutil
