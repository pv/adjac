!! NOTE: this file is autogenerated from adjac_fft.f95.in: do not edit manually
module adjac_fft
  use adjac

  private
  public fft, ifft

  interface fft
     module procedure fft_a, fft_d, fft_b, fft_z
  end interface fft

  interface ifft
     module procedure ifft_a, ifft_d, ifft_b, ifft_z
  end interface ifft

  type(adjac_double), dimension(:), allocatable :: ch_a
  double precision, dimension(:), allocatable :: wa, ch_d
  integer, dimension(15) :: ifac
contains
  subroutine fft_a(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    type(adjac_double), dimension(:), contiguous, intent(inout) :: q
    integer :: n

    n = size(q)/2
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_a) .or. size(ch_a) .ne. 2*n) then
       if (allocated(ch_a)) deallocate(ch_a)
       allocate(ch_a(2*n))
    end if
    call zfftf1a(n, q, ch_a, wa, ifac)
  end subroutine fft_a

  subroutine ifft_a(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    type(adjac_double), dimension(:), contiguous, intent(inout) :: q
    integer :: i, n

    n = size(q)/2
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_a) .or. size(ch_a) .ne. 2*n) then
       if (allocated(ch_a)) deallocate(ch_a)
       allocate(ch_a(2*n))
    end if
    call zfftb1a(n, q, ch_a, wa, ifac)

    do i = 1, size(q)
       q(i) = q(i) / n
    end do
  end subroutine ifft_a
  subroutine fft_d(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    double precision, dimension(:), contiguous, intent(inout) :: q
    integer :: n

    n = size(q)/2
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_d) .or. size(ch_d) .ne. 2*n) then
       if (allocated(ch_d)) deallocate(ch_d)
       allocate(ch_d(2*n))
    end if
    call zfftf1d(n, q, ch_d, wa, ifac)
  end subroutine fft_d

  subroutine ifft_d(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    double precision, dimension(:), contiguous, intent(inout) :: q
    integer :: i, n

    n = size(q)/2
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_d) .or. size(ch_d) .ne. 2*n) then
       if (allocated(ch_d)) deallocate(ch_d)
       allocate(ch_d(2*n))
    end if
    call zfftb1d(n, q, ch_d, wa, ifac)

    do i = 1, size(q)
       q(i) = q(i) / n
    end do
  end subroutine ifft_d
  subroutine fft_b(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    type(adjac_complex), dimension(:), contiguous, intent(inout) :: q
    integer :: n

    n = size(q)
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_a) .or. size(ch_a) .ne. 2*n) then
       if (allocated(ch_a)) deallocate(ch_a)
       allocate(ch_a(2*n))
    end if
    call zfftf1a(n, q, ch_a, wa, ifac)
  end subroutine fft_b

  subroutine ifft_b(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    type(adjac_complex), dimension(:), contiguous, intent(inout) :: q
    integer :: i, n

    n = size(q)
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_a) .or. size(ch_a) .ne. 2*n) then
       if (allocated(ch_a)) deallocate(ch_a)
       allocate(ch_a(2*n))
    end if
    call zfftb1a(n, q, ch_a, wa, ifac)

    do i = 1, size(q)
       q(i) = q(i) / n
    end do
  end subroutine ifft_b
  subroutine fft_z(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    double complex, dimension(:), contiguous, intent(inout) :: q
    integer :: n

    n = size(q)
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_d) .or. size(ch_d) .ne. 2*n) then
       if (allocated(ch_d)) deallocate(ch_d)
       allocate(ch_d(2*n))
    end if
    call zfftf1d(n, q, ch_d, wa, ifac)
  end subroutine fft_z

  subroutine ifft_z(q)
    use adjac
    implicit none
    external :: zffti1, zfftf1
    double complex, dimension(:), contiguous, intent(inout) :: q
    integer :: i, n

    n = size(q)
    if (n.le.1) return

    if (.not.allocated(wa) .or. size(wa) .ne. 2*n) then
       if (allocated(wa)) deallocate(wa)
       allocate(wa(2*n))
       call zffti1(n, wa, ifac)
    end if
    if (.not.allocated(ch_d) .or. size(ch_d) .ne. 2*n) then
       if (allocated(ch_d)) deallocate(ch_d)
       allocate(ch_d(2*n))
    end if
    call zfftb1d(n, q, ch_d, wa, ifac)

    do i = 1, size(q)
       q(i) = q(i) / n
    end do
  end subroutine ifft_z
end module adjac_fft