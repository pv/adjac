! -*-f90-*-
module testutil

  interface assert_equal
     module procedure assert_equal_i, assert_equal_d, assert_equal_z
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

end module testutil
