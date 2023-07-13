program test_naca_geometry
  use kinds
  use constants, only: pi
  use naca_geometry
  use reference_geometry
  use tolerances
  use test_support
  implicit none

  real(real64), parameter :: stations(6) = [0.0125_real64, 0.05_real64, 0.2_real64, 0.3_real64, 0.6_real64, 1.0_real64]
  real(real64) :: x(256), y(256), err, face, expected, a_l0, th, slope, dth, mm, pp
  integer :: ierr, i, ni, nq

  ! independent thickness evaluation
  err = 0.0_real64
  do i = 1, 6
    err = max(err, abs(thickness('0012', stations(i)) - ref_thickness(0.12_real64, stations(i))))
  end do
  call check_close(err, 0.0_real64, 1.0e-12_real64, 'thickness vs independent evaluator')
  call check_close(thickness('0012', 0.3_real64), 0.0600_real64, 1.0e-4_real64, 'NACA 0012 max half thickness')

  ! NACA 0012 ring: clockwise, mirror symmetric, node 1 on the lower side
  ni = 256
  call naca_ring('0012', ni, x, y, ierr)
  call check(ierr == ierr_ok, 'naca_ring ok')
  call check(shoelace(x(1:ni), y(1:ni)) < 0.0_real64, 'ring is clockwise')
  call check(y(1) < 0.0_real64 .and. y(ni) > 0.0_real64, 'node 1 lower, node ni upper')
  err = 0.0_real64
  do i = 1, ni
    err = max(err, abs(x(i) - x(ni + 1 - i)), abs(y(i) + y(ni + 1 - i)))
  end do
  call check_close(err, 0.0_real64, 1.0e-14_real64, 'symmetric section mirror symmetry')
  face = abs(y(ni) - y(1))
  expected = 2.0_real64 * 0.6_real64 * (0.2969_real64 - 0.1260_real64 - 0.3516_real64 + 0.2843_real64 - 0.1015_real64)
  call check_close(face / expected, 1.0_real64, tol_naca_face_rel, 'blunt face thickness (relative)')
  call check_close(x(1), 1.0_real64, 1.0e-4_real64, 'trailing edge near x=1')

  ! cambered 2412: upper x differs from lower x; camber line midway
  call naca_ring('2412', ni, x, y, ierr)
  call check(ierr == ierr_ok, 'naca_ring 2412 ok')
  call check(abs(x(ni) - x(1)) > 0.0_real64 .or. abs(y(ni) + y(1)) > 0.0_real64, '2412 is not mirror symmetric')

  ! zero-lift angle of 2412 from camber slope: -(1/pi) int_0^pi slope (cos t - 1) dt  (midpoint rule)
  mm = 0.02_real64; pp = 0.4_real64
  nq = 20000
  a_l0 = 0.0_real64
  dth = pi / nq
  do i = 1, nq
    th = (i - 0.5_real64) * dth
    slope = ref_camber_slope(mm, pp, 0.5_real64 * (1.0_real64 - cos(th)))
    a_l0 = a_l0 - slope * (cos(th) - 1.0_real64) * dth / pi
  end do
  print '(a,f8.4,a)', '  2412 thin-airfoil zero-lift angle ', a_l0 * 180.0_real64 / pi, ' deg'
  call check_close(a_l0 * 180.0_real64 / pi, -2.077_real64, 0.01_real64, '2412 alpha_L0 near -2.08 deg')
  call check_close(ref_camber(mm, pp, 0.4_real64), camber('2412', 0.4_real64), 1.0e-14_real64, 'camber vs reference')

  ! error cases
  call naca_ring('12', 64, x, y, ierr);   call check(ierr == ierr_bad_input, 'short code rejected')
  call naca_ring('ab12', 64, x, y, ierr); call check(ierr == ierr_bad_input, 'non-digit rejected')
  call naca_ring('2400', 64, x, y, ierr); call check(ierr == ierr_bad_input, 'zero thickness rejected')
  call naca_ring('2012', 64, x, y, ierr); call check(ierr == ierr_bad_input, 'camber without position rejected')
  call naca_ring('0012', 63, x, y, ierr); call check(ierr == ierr_bad_input, 'odd ni rejected')
  call naca_ring('0012', 16, x, y, ierr); call check(ierr == ierr_bad_input, 'small ni rejected')

  call finish()
end program test_naca_geometry
