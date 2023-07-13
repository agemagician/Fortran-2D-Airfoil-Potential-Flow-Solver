! Test-only reference shapes and exact answers, independent of src/ geometry code.
module reference_geometry
  use kinds, only: real64
  use constants, only: pi, center_x, center_y
  implicit none
  public
  real(real64), parameter :: jk_a = 0.25_real64, jk_cx = -0.025_real64, jk_cy = 0.02_real64
  real(real64), parameter :: jk_rc = sqrt((jk_a - jk_cx)**2 + jk_cy**2)
contains

  ! Independent NACA thickness in Horner form, t = thickness fraction.
  function ref_thickness(t, xc) result(yt)
    real(real64), intent(in) :: t, xc
    real(real64) :: yt
    yt = 5.0_real64 * t * (0.2969_real64 * sqrt(xc) + xc * (-0.1260_real64 + xc * (-0.3516_real64 &
         + xc * (0.2843_real64 + xc * (-0.1015_real64)))))
  end function ref_thickness

  ! Independent NACA camber line (m, p as fractions) and slope.
  function ref_camber(m, p, xc) result(yc)
    real(real64), intent(in) :: m, p, xc
    real(real64) :: yc
    if (m == 0.0_real64) then
      yc = 0.0_real64
    else if (xc < p) then
      yc = m / p**2 * (2.0_real64 * p * xc - xc**2)
    else
      yc = m / (1.0_real64 - p)**2 * (1.0_real64 - 2.0_real64 * p + 2.0_real64 * p * xc - xc**2)
    end if
  end function ref_camber

  function ref_camber_slope(m, p, xc) result(s)
    real(real64), intent(in) :: m, p, xc
    real(real64) :: s
    if (m == 0.0_real64) then
      s = 0.0_real64
    else if (xc < p) then
      s = 2.0_real64 * m / p**2 * (p - xc)
    else
      s = 2.0_real64 * m / (1.0_real64 - p)**2 * (p - xc)
    end if
  end function ref_camber_slope

  ! Clockwise circle of radius a about the outer-circle centre, first node at angle -pi/ni.
  subroutine circle_ring(a, ni, x, y)
    real(real64), intent(in) :: a
    integer, intent(in) :: ni
    real(real64), intent(out) :: x(ni), y(ni)
    integer :: i
    real(real64) :: ang
    do i = 1, ni
      ang = -pi / ni - 2.0_real64 * pi * (i - 1) / ni
      x(i) = center_x + a * cos(ang)
      y(i) = center_y + a * sin(ang)
    end do
  end subroutine circle_ring

  ! Signed shoelace area (negative for a clockwise ring).
  function shoelace(x, y) result(s)
    real(real64), intent(in) :: x(:), y(:)
    real(real64) :: s
    integer :: i, j, n
    n = size(x)
    s = 0.0_real64
    do i = 1, n
      j = mod(i, n) + 1
      s = s + x(i) * y(j) - x(j) * y(i)
    end do
    s = 0.5_real64 * s
  end function shoelace

  ! ---- Joukowski airfoil: z = zeta + a^2/zeta, circle through zeta = a, centre (-0.025, 0.02) ----
  function jk_xy(th, xmin, sc) result(p)
    real(real64), intent(in) :: th, xmin, sc
    real(real64) :: p(2)
    complex(real64) :: zeta, z
    zeta = cmplx(jk_cx + jk_rc * cos(th), jk_cy + jk_rc * sin(th), real64)
    z = zeta + jk_a**2 / zeta
    p = [(real(z) - xmin) * sc, aimag(z) * sc]
  end function jk_xy

  function jk_x(th, xmin, sc) result(xq)
    real(real64), intent(in) :: th, xmin, sc
    real(real64) :: xq, p(2)
    p = jk_xy(th, xmin, sc)
    xq = p(1)
  end function jk_x

  ! Leading-edge parameter (minimum x), scale to unit x-extent, and trailing-edge parameter.
  subroutine jk_setup(th_te, th_le, xmin, sc)
    real(real64), intent(out) :: th_te, th_le, xmin, sc
    real(real64) :: lo, hi, m1, m2, gr
    integer :: it
    th_te = atan2(-jk_cy, jk_a - jk_cx)
    lo = th_te - 2.0_real64 * pi
    hi = th_te
    gr = 0.5_real64 * (sqrt(5.0_real64) - 1.0_real64)
    do it = 1, 200
      m1 = hi - gr * (hi - lo)
      m2 = lo + gr * (hi - lo)
      if (jk_x(m1, 0.0_real64, 1.0_real64) < jk_x(m2, 0.0_real64, 1.0_real64)) then
        hi = m2
      else
        lo = m1
      end if
    end do
    th_le = 0.5_real64 * (lo + hi)
    xmin = jk_x(th_le, 0.0_real64, 1.0_real64)
    sc = 1.0_real64 / (jk_x(th_te, 0.0_real64, 1.0_real64) - xmin)
  end subroutine jk_setup

  ! y on the lower or upper branch at scaled abscissa xq.
  function jk_y(xq, lower, th_te, th_le, xmin, sc) result(yq)
    real(real64), intent(in) :: xq, th_te, th_le, xmin, sc
    logical, intent(in) :: lower
    real(real64) :: yq, lo, hi, mid, xl, p(2)
    integer :: it
    if (lower) then
      lo = th_le; hi = th_te        ! x increases with theta on the lower branch
    else
      lo = th_te - 2.0_real64 * pi; hi = th_le   ! x decreases with theta on the upper branch
    end if
    do it = 1, 100
      mid = 0.5_real64 * (lo + hi)
      xl = jk_x(mid, xmin, sc)
      if ((xl < xq) .eqv. lower) then
        lo = mid
      else
        hi = mid
      end if
    end do
    p = jk_xy(0.5_real64 * (lo + hi), xmin, sc)
    yq = p(2)
  end function jk_y

  ! Blunt-face Joukowski ring (clockwise, node 1 = lower corner), face cut where the vertical gap is e.
  subroutine joukowski_ring(e, ni, x, y)
    real(real64), intent(in) :: e
    integer, intent(in) :: ni
    real(real64), intent(out) :: x(ni), y(ni)
    real(real64) :: th_te, th_le, xmin, sc, lo, hi, xc, phi, xk
    integer :: m, k, it
    call jk_setup(th_te, th_le, xmin, sc)
    lo = 0.5_real64; hi = 1.0_real64
    do it = 1, 100
      xc = 0.5_real64 * (lo + hi)
      if (jk_y(xc, .false., th_te, th_le, xmin, sc) - jk_y(xc, .true., th_te, th_le, xmin, sc) > e) then
        lo = xc
      else
        hi = xc
      end if
    end do
    xc = 0.5_real64 * (lo + hi)
    m = ni / 2
    do k = 1, m
      phi = pi * (k - 0.5_real64) / (m - 0.5_real64)
      xk = xc * 0.5_real64 * (1.0_real64 - cos(phi))
      x(m + k) = xk
      y(m + k) = jk_y(xk, .false., th_te, th_le, xmin, sc)
      x(m - k + 1) = xk
      y(m - k + 1) = jk_y(xk, .true., th_te, th_le, xmin, sc)
    end do
  end subroutine joukowski_ring

  ! Exact CL of the sharp-edged (unit x-extent) Joukowski airfoil: 2 Gamma / U, Gamma = 4 pi U Rc sc sin(alpha + beta).
  function joukowski_cl(alpha_rad) result(cl)
    real(real64), intent(in) :: alpha_rad
    real(real64) :: cl, th_te, th_le, xmin, sc, beta
    call jk_setup(th_te, th_le, xmin, sc)
    beta = asin(jk_cy / jk_rc)
    cl = 2.0_real64 * 4.0_real64 * pi * jk_rc * sc * sin(alpha_rad + beta)
  end function joukowski_cl

end module reference_geometry
