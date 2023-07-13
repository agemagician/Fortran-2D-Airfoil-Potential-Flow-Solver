! O-grid builder, metrics and quality gate. Works on any closed clockwise ring.
module grid
  use kinds, only: real64, ierr_ok, ierr_bad_input, ierr_spacing, ierr_quality
  use constants, only: pi, center_x, center_y, chord, max_skew_deg, max_aspect, min_area_ratio, te_wedge_rad
  implicit none
  private
  public :: grid_record, metrics_record, quality_report
  public :: build_ogrid, grid_metrics, check_quality, cell_area, quality_failure

  type :: grid_record
    integer :: ni = 0, nj = 0
    real(real64), allocatable :: x(:, :), y(:, :)
  end type grid_record

  type :: metrics_record
    real(real64), allocatable, dimension(:, :) :: x_xi, x_eta, y_xi, y_eta, alpha_g, beta_g, gamma_g, jac
  end type metrics_record

  type :: quality_report
    real(real64) :: min_area = 0.0_real64, max_skew_deg = 0.0_real64, max_aspect = 0.0_real64
    real(real64) :: area_ratio = 0.0_real64, min_jac = 0.0_real64
    integer :: first_bad_i = 0, first_bad_j = 0
    character(len=16) :: limit = ''   ! name of the first failed limit, blank if none
  end type quality_report

contains

  ! exp(x) - 1 without cancellation for small x (expm1 is not in Fortran 2008).
  pure function em1(x) result(e)
    real(real64), intent(in) :: x
    real(real64) :: e
    if (abs(x) < 1.0e-5_real64) then
      e = x * (1.0_real64 + 0.5_real64 * x * (1.0_real64 + x / 3.0_real64))
    else
      e = exp(x) - 1.0_real64
    end if
  end function em1

  ! g(r) = h ((1+r)^n - 1)/r : distance to the outermost node for growth ratio q = 1 + r.
  pure function span(h, n, r) result(g)
    real(real64), intent(in) :: h, r
    integer, intent(in) :: n
    real(real64) :: g
    if (r < 1.0e-12_real64) then
      g = h * n
    else
      g = h * em1(n * log(1.0_real64 + r)) / r
    end if
  end function span

  ! Ring (x, y) clockwise with ni even >= 32; nj >= 16. Straight radial lines, geometric spacing.
  subroutine build_ogrid(x, y, nj, router, first_cell, g, ierr)
    real(real64), intent(in) :: x(:), y(:), router, first_cell
    integer, intent(in) :: nj
    type(grid_record), intent(out) :: g
    integer, intent(out) :: ierr
    real(real64), allocatable :: s(:)
    real(real64) :: seg, total, f, theta0, turn, phi, xo, yo, len_i, h, lo, hi, r, t, d
    integer :: ni, i, j, it

    ierr = ierr_ok
    ni = size(x)
    if (ni < 32 .or. mod(ni, 2) /= 0 .or. size(y) /= ni .or. nj < 16 .or. router <= 0.0_real64 &
        .or. first_cell <= 0.0_real64) then
      ierr = ierr_bad_input
      return
    end if
    h = first_cell * chord

    allocate (s(ni))
    s(1) = 0.0_real64
    do i = 2, ni
      s(i) = s(i - 1) + hypot(x(i) - x(i - 1), y(i) - y(i - 1))
    end do
    seg = hypot(x(1) - x(ni), y(1) - y(ni))
    total = s(ni) + seg
    f = seg / total
    ! Half-angle of the outer-ring gap behind the body. A sharp turn at the face corner (airfoil) needs a wide
    ! wedge so the first radial lines leave the corner inside the exterior angle; smooth rings keep pi*f.
    turn = abs(atan2((x(1) - x(ni)) * (y(2) - y(1)) - (y(1) - y(ni)) * (x(2) - x(1)), &
                     (x(1) - x(ni)) * (x(2) - x(1)) + (y(1) - y(ni)) * (y(2) - y(1))))
    theta0 = pi * f
    if (turn > 0.25_real64 * pi) theta0 = te_wedge_rad
    s = s / s(ni)    ! s(ni) = 1: body nodes span the whole outer arc 2 pi - 2 theta0

    g%ni = ni
    g%nj = nj
    allocate (g%x(ni, nj), g%y(ni, nj))
    do i = 1, ni
      phi = -theta0 - (2.0_real64 * pi - 2.0_real64 * theta0) * s(i)
      xo = center_x + router * chord * cos(phi)
      yo = center_y + router * chord * sin(phi)
      len_i = hypot(xo - x(i), yo - y(i))
      if (h * (nj - 1) > len_i * (1.0_real64 + 1.0e-12_real64)) then
        ierr = ierr_spacing
        write (*, '(a,es10.3,a,es10.3)') 'first_cell too large for nj and router; largest allowed first_cell = ', &
          len_i / (nj - 1) / chord, ' at ring node; given ', first_cell
        return
      end if
      lo = 0.0_real64
      hi = 1.0_real64
      do while (span(h, nj - 1, hi) < len_i)
        hi = 2.0_real64 * hi
      end do
      do it = 1, 200
        r = 0.5_real64 * (lo + hi)
        if (span(h, nj - 1, r) < len_i) then
          lo = r
        else
          hi = r
        end if
      end do
      r = 0.5_real64 * (lo + hi)
      do j = 1, nj
        if (j == nj) then
          t = 1.0_real64
        else if (r < 1.0e-12_real64) then
          t = h * (j - 1) / len_i
        else
          d = h * em1((j - 1) * log(1.0_real64 + r)) / r
          t = d / len_i
        end if
        g%x(i, j) = x(i) + t * (xo - x(i))
        g%y(i, j) = y(i) + t * (yo - y(i))
      end do
    end do
  end subroutine build_ogrid

  subroutine grid_metrics(g, met)
    type(grid_record), intent(in) :: g
    type(metrics_record), intent(out) :: met
    integer :: ni, nj, i, j, ip, im
    ni = g%ni
    nj = g%nj
    allocate (met%x_xi(ni, nj), met%x_eta(ni, nj), met%y_xi(ni, nj), met%y_eta(ni, nj))
    allocate (met%alpha_g(ni, nj), met%beta_g(ni, nj), met%gamma_g(ni, nj), met%jac(ni, nj))
    do j = 1, nj
      do i = 1, ni
        ip = mod(i, ni) + 1
        im = mod(i - 2 + ni, ni) + 1
        met%x_xi(i, j) = 0.5_real64 * (g%x(ip, j) - g%x(im, j))
        met%y_xi(i, j) = 0.5_real64 * (g%y(ip, j) - g%y(im, j))
        if (j == 1) then
          met%x_eta(i, j) = 0.5_real64 * (-3.0_real64 * g%x(i, 1) + 4.0_real64 * g%x(i, 2) - g%x(i, 3))
          met%y_eta(i, j) = 0.5_real64 * (-3.0_real64 * g%y(i, 1) + 4.0_real64 * g%y(i, 2) - g%y(i, 3))
        else if (j == nj) then
          met%x_eta(i, j) = 0.5_real64 * (3.0_real64 * g%x(i, nj) - 4.0_real64 * g%x(i, nj - 1) + g%x(i, nj - 2))
          met%y_eta(i, j) = 0.5_real64 * (3.0_real64 * g%y(i, nj) - 4.0_real64 * g%y(i, nj - 1) + g%y(i, nj - 2))
        else
          met%x_eta(i, j) = 0.5_real64 * (g%x(i, j + 1) - g%x(i, j - 1))
          met%y_eta(i, j) = 0.5_real64 * (g%y(i, j + 1) - g%y(i, j - 1))
        end if
      end do
    end do
    met%alpha_g = met%x_eta**2 + met%y_eta**2
    met%beta_g = met%x_xi * met%x_eta + met%y_xi * met%y_eta
    met%gamma_g = met%x_xi**2 + met%y_xi**2
    met%jac = met%x_xi * met%y_eta - met%x_eta * met%y_xi
  end subroutine grid_metrics

  ! Cell (i,j) has corners (i,j), (i+1,j), (i+1,j+1), (i,j+1); i wraps.
  subroutine cell_corners(g, i, j, cx, cy)
    type(grid_record), intent(in) :: g
    integer, intent(in) :: i, j
    real(real64), intent(out) :: cx(4), cy(4)
    integer :: ip
    ip = mod(i, g%ni) + 1
    cx = [g%x(i, j), g%x(ip, j), g%x(ip, j + 1), g%x(i, j + 1)]
    cy = [g%y(i, j), g%y(ip, j), g%y(ip, j + 1), g%y(i, j + 1)]
  end subroutine cell_corners

  function cell_area(g, i, j) result(a)
    type(grid_record), intent(in) :: g
    integer, intent(in) :: i, j
    real(real64) :: a, cx(4), cy(4)
    call cell_corners(g, i, j, cx, cy)
    a = 0.5_real64 * ((cx(1) - cx(3)) * (cy(2) - cy(4)) - (cx(2) - cx(4)) * (cy(1) - cy(3)))
  end function cell_area

  subroutine cell_shape(g, i, j, skew, aspect)
    type(grid_record), intent(in) :: g
    integer, intent(in) :: i, j
    real(real64), intent(out) :: skew, aspect
    real(real64) :: cx(4), cy(4), ax, ay, bx, by, na, nb, c, ang, e(4), lx, ly
    integer :: k, kp, km
    call cell_corners(g, i, j, cx, cy)
    skew = 0.0_real64
    do k = 1, 4
      kp = mod(k, 4) + 1
      km = mod(k + 2, 4) + 1
      ax = cx(kp) - cx(k); ay = cy(kp) - cy(k)
      bx = cx(km) - cx(k); by = cy(km) - cy(k)
      na = hypot(ax, ay); nb = hypot(bx, by)
      if (na == 0.0_real64 .or. nb == 0.0_real64) then
        skew = 90.0_real64
      else
        c = max(-1.0_real64, min(1.0_real64, (ax * bx + ay * by) / (na * nb)))
        ang = acos(c) * 180.0_real64 / pi
        skew = max(skew, abs(90.0_real64 - ang))
      end if
    end do
    do k = 1, 4
      kp = mod(k, 4) + 1
      e(k) = hypot(cx(kp) - cx(k), cy(kp) - cy(k))
    end do
    lx = 0.5_real64 * (e(1) + e(3))
    ly = 0.5_real64 * (e(2) + e(4))
    if (min(lx, ly) <= 0.0_real64) then
      aspect = huge(1.0_real64)
    else
      aspect = max(lx, ly) / min(lx, ly)
    end if
  end subroutine cell_shape

  pure function quality_failure(rep) result(name)
    type(quality_report), intent(in) :: rep
    character(len=16) :: name
    name = ''
    if (.not. (rep%min_area > 0.0_real64)) then
      name = 'min_area'
    else if (.not. (rep%min_jac > 0.0_real64)) then
      name = 'min_jacobian'
    else if (.not. (rep%max_skew_deg <= max_skew_deg)) then
      name = 'max_skew'
    else if (.not. (rep%max_aspect <= max_aspect)) then
      name = 'max_aspect'
    else if (.not. (rep%area_ratio >= min_area_ratio)) then
      name = 'area_ratio'
    end if
  end function quality_failure

  recursive subroutine qsort(a)
    real(real64), intent(inout) :: a(:)
    real(real64) :: p, tmp
    integer :: i, j
    if (size(a) < 2) return
    p = a((size(a) + 1) / 2)
    i = 1
    j = size(a)
    do
      do while (a(i) < p)
        i = i + 1
      end do
      do while (a(j) > p)
        j = j - 1
      end do
      if (i >= j) exit
      tmp = a(i); a(i) = a(j); a(j) = tmp
      i = i + 1
      j = j - 1
    end do
    call qsort(a(1:j))
    call qsort(a(j + 1:))
  end subroutine qsort

  subroutine check_quality(g, rep, ierr)
    type(grid_record), intent(in) :: g
    type(quality_report), intent(out) :: rep
    integer, intent(out) :: ierr
    type(metrics_record) :: met
    real(real64), allocatable :: areas(:)
    real(real64) :: sk, asp, a, median
    integer :: i, j, n, nc

    ierr = ierr_ok
    call grid_metrics(g, met)
    nc = g%ni * (g%nj - 1)
    allocate (areas(nc))
    rep%max_skew_deg = 0.0_real64
    rep%max_aspect = 0.0_real64
    n = 0
    do j = 1, g%nj - 1
      do i = 1, g%ni
        n = n + 1
        areas(n) = cell_area(g, i, j)
        call cell_shape(g, i, j, sk, asp)
        rep%max_skew_deg = max(rep%max_skew_deg, sk)
        rep%max_aspect = max(rep%max_aspect, asp)
      end do
    end do
    rep%min_area = minval(areas)
    rep%min_jac = minval(met%jac)
    call qsort(areas)
    median = areas((nc + 1) / 2)
    rep%area_ratio = rep%min_area / median
    rep%limit = quality_failure(rep)
    if (len_trim(rep%limit) == 0) return

    ierr = ierr_quality
    ! first offending cell (lowest j, then lowest i) for the failed limit
    scan: do j = 1, g%nj - 1
      do i = 1, g%ni
        a = cell_area(g, i, j)
        call cell_shape(g, i, j, sk, asp)
        select case (trim(rep%limit))
        case ('min_area')
          if (.not. (a > 0.0_real64)) exit scan
        case ('min_jacobian')
          if (.not. (met%jac(i, j) > 0.0_real64)) exit scan
        case ('max_skew')
          if (.not. (sk <= max_skew_deg)) exit scan
        case ('max_aspect')
          if (.not. (asp <= max_aspect)) exit scan
        case default
          if (.not. (a / median >= min_area_ratio)) exit scan
        end select
      end do
    end do scan
    rep%first_bad_i = min(i, g%ni)
    rep%first_bad_j = min(j, g%nj - 1)
  end subroutine check_quality

end module grid
