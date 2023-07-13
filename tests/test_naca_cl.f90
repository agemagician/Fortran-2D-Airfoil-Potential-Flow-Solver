program test_naca_cl
  use kinds
  use constants, only: pi, cl_gap_flag
  use naca_geometry
  use grid
  use potential_flow
  use postproc
  use reference_geometry
  use tolerances
  use test_support
  implicit none

  integer, parameter :: ni = 256, nj = 96
  real(real64) :: cl_c, cl_p, gap, cl_m, cl_p2, a_l0, th, dth, thin, cl_thick, slope
  integer :: ierr, i
  type(force_report) :: rep

  call run('0012', 0.0_real64, cl_c, cl_p, gap, ierr)
  call check(ierr == ierr_ok, 'NACA 0012 alpha=0 solve ok')
  call check_close(cl_p, 0.0_real64, tol_zero_cl, 'NACA 0012 alpha=0 CL from pressure')

  call run('0012', 5.0_real64, cl_c, cl_p, gap, ierr)
  cl_thick = 2.0_real64 * pi * (1.0_real64 + 0.77_real64 * 0.12_real64) * 5.0_real64 * pi / 180.0_real64
  print '(a,f8.4,a,f8.4,a,f8.4)', '  0012 at 5 deg: CL circ ', cl_c, '  CL press ', cl_p, '  rule-of-thumb ', cl_thick
  if (abs(cl_c / cl_thick - 1.0_real64) > 0.05_real64) print '(a)', '  WARNING: CL outside 5% of the thickness-corrected thin-airfoil value (not a gate)'

  ! 2412 zero-lift angle: two solves, linear interpolation, against thin-airfoil integral
  call run('2412', -2.0_real64, cl_m, cl_p, gap, ierr)
  call run('2412', 2.0_real64, cl_p2, cl_p, gap, ierr)
  a_l0 = -2.0_real64 + 4.0_real64 * (0.0_real64 - cl_m) / (cl_p2 - cl_m)
  thin = 0.0_real64
  dth = pi / 20000
  do i = 1, 20000
    th = (i - 0.5_real64) * dth
    slope = ref_camber_slope(0.02_real64, 0.4_real64, 0.5_real64 * (1.0_real64 - cos(th)))
    thin = thin - slope * (cos(th) - 1.0_real64) * dth / pi
  end do
  print '(a,f8.4,a,f8.4)', '  2412 zero-lift angle: solver ', a_l0, ' deg, thin-airfoil ', thin * 180.0_real64 / pi
  call check_close(a_l0, thin * 180.0_real64 / pi, tol_alpha_l0_deg, '2412 zero-lift angle vs thin-airfoil')

  call run('2412', 5.0_real64, cl_c, cl_p, gap, ierr)
  call pressure_forces_print(gap)

  call finish()

contains

  subroutine run(code, alpha_deg, cl_circ, cl_press, cl_gap, ierr_out)
    character(*), intent(in) :: code
    real(real64), intent(in) :: alpha_deg
    real(real64), intent(out) :: cl_circ, cl_press, cl_gap
    integer, intent(out) :: ierr_out
    real(real64) :: x(ni), y(ni), cp(ni)
    real(real64), allocatable :: u(:, :), v(:, :)
    type(grid_record) :: g
    type(metrics_record) :: met
    type(solver_params) :: par
    type(flow_result) :: r

    call naca_ring(code, ni, x, y, ierr_out)
    call build_ogrid(x, y, nj, 30.0_real64, 1.0e-3_real64, g, ierr_out)
    call grid_metrics(g, met)
    allocate (u(ni, nj), v(ni, nj))
    par%alpha_rad = alpha_deg * pi / 180.0_real64
    par%maxiter = 5000
    par%tol = 1.0e-13_real64
    call solve_flow(g, met, par, r)
    ierr_out = r%ierr
    cl_circ = 2.0_real64 * r%gamma
    call velocity_field(met, r%psi, u, v)
    call surface_cp(u, v, cp)
    call pressure_forces(g, cp, par%alpha_rad, cl_circ, rep)
    cl_press = rep%cl_press
    cl_gap = rep%cl_gap_rel
    print '(a,a,f6.1,a,f9.5,a,f9.5,a,es10.3,a,es10.3)', '  ', code, alpha_deg, ' deg  CL circ ', cl_circ, '  CL press ', &
      cl_press, '  gap ', cl_gap, '  CD ', rep%cd
  end subroutine run

  subroutine pressure_forces_print(gap_rel)
    real(real64), intent(in) :: gap_rel
    if (gap_rel > cl_gap_flag) print '(a,es10.3,a,es10.3)', 'FLAG: cl gap ', gap_rel, ' exceeds ', cl_gap_flag
  end subroutine pressure_forces_print

end program test_naca_cl
