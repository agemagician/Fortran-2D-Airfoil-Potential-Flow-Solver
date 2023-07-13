! NACA 4-digit surface ring with a blunt trailing edge (-0.1015 thickness coefficient).
module naca_geometry
  use kinds, only: real64, ierr_ok, ierr_bad_input
  use constants, only: pi
  implicit none
  private
  public :: naca_ring, thickness, camber

contains

  ! Parses a 4-digit code. ok is false for non-digits, zero thickness, or camber with zero position.
  subroutine parse_code(code, mc, pc, tc, ok)
    character(*), intent(in) :: code
    real(real64), intent(out) :: mc, pc, tc
    logical, intent(out) :: ok
    integer :: d(4), i
    mc = 0.0_real64; pc = 0.0_real64; tc = 0.0_real64
    ok = len_trim(code) == 4
    if (.not. ok) return
    do i = 1, 4
      d(i) = iachar(code(i:i)) - iachar('0')
      if (d(i) < 0 .or. d(i) > 9) then
        ok = .false.
        return
      end if
    end do
    mc = d(1) / 100.0_real64
    pc = d(2) / 10.0_real64
    tc = (10 * d(3) + d(4)) / 100.0_real64
    ok = tc > 0.0_real64 .and. (d(1) == 0 .or. d(2) > 0)
  end subroutine parse_code

  pure function thickness(code, xc) result(yt)
    character(*), intent(in) :: code
    real(real64), intent(in) :: xc
    real(real64) :: yt, tc
    tc = (10 * (iachar(code(3:3)) - iachar('0')) + iachar(code(4:4)) - iachar('0')) / 100.0_real64
    yt = 5.0_real64 * tc * (0.2969_real64 * sqrt(xc) - 0.1260_real64 * xc - 0.3516_real64 * xc**2 &
                            + 0.2843_real64 * xc**3 - 0.1015_real64 * xc**4)
  end function thickness

  pure function camber(code, xc) result(yc)
    character(*), intent(in) :: code
    real(real64), intent(in) :: xc
    real(real64) :: yc, mc, pc
    mc = (iachar(code(1:1)) - iachar('0')) / 100.0_real64
    pc = (iachar(code(2:2)) - iachar('0')) / 10.0_real64
    if (mc == 0.0_real64) then
      yc = 0.0_real64
    else if (xc < pc) then
      yc = mc / pc**2 * (2.0_real64 * pc * xc - xc**2)
    else
      yc = mc / (1.0_real64 - pc)**2 * ((1.0_real64 - 2.0_real64 * pc) + 2.0_real64 * pc * xc - xc**2)
    end if
  end function camber

  pure function camber_slope(code, xc) result(dyc)
    character(*), intent(in) :: code
    real(real64), intent(in) :: xc
    real(real64) :: dyc, mc, pc
    mc = (iachar(code(1:1)) - iachar('0')) / 100.0_real64
    pc = (iachar(code(2:2)) - iachar('0')) / 10.0_real64
    if (mc == 0.0_real64) then
      dyc = 0.0_real64
    else if (xc < pc) then
      dyc = 2.0_real64 * mc / pc**2 * (pc - xc)
    else
      dyc = 2.0_real64 * mc / (1.0_real64 - pc)**2 * (pc - xc)
    end if
  end function camber_slope

  ! Clockwise ring: node 1 = lower trailing-edge corner, node ni = upper corner.
  subroutine naca_ring(code, ni, x, y, ierr)
    character(*), intent(in) :: code
    integer, intent(in) :: ni
    real(real64), intent(out) :: x(:), y(:)
    integer, intent(out) :: ierr
    real(real64) :: mc, pc, tc, phi, xs, yt, yc, th
    integer :: m, k
    logical :: ok

    ierr = ierr_ok
    x = 0.0_real64; y = 0.0_real64
    call parse_code(code, mc, pc, tc, ok)
    if (.not. ok .or. ni < 32 .or. mod(ni, 2) /= 0 .or. size(x) < ni .or. size(y) < ni) then
      ierr = ierr_bad_input
      return
    end if
    m = ni / 2
    do k = 1, m
      phi = pi * (k - 0.5_real64) / (m - 0.5_real64)
      xs = 0.5_real64 * (1.0_real64 - cos(phi))
      yt = thickness(code, xs)
      yc = camber(code, xs)
      th = atan(camber_slope(code, xs))
      x(m + k) = xs - yt * sin(th)
      y(m + k) = yc + yt * cos(th)
      x(m - k + 1) = xs + yt * sin(th)
      y(m - k + 1) = yc - yt * cos(th)
    end do
  end subroutine naca_ring

end module naca_geometry
