program test_cylinder
  use kinds
  use constants, only: pi, uinf, center_x, center_y
  use grid
  use potential_flow
  use postproc
  use reference_geometry
  use tolerances
  use test_support
  implicit none

  real(real64), parameter :: a_cyl = 0.5_real64, router = 30.0_real64, gamma_target = 1.0_real64
  integer, parameter :: ni = 256, nj = 96
  type(grid_record) :: g
  type(metrics_record) :: met
  type(solver_params) :: par
  type(flow_result) :: r0, r1
  real(real64), allocatable :: rx(:), ry(:), u(:, :), v(:, :), cp(:), cp_exact(:), q(:), qa(:)
  real(real64) :: gam_aux, big_r, k, gamma_aux_exact, th, speed_ex, c, gam, top0, top1
  integer :: ierr, i, itop

  allocate (rx(ni), ry(ni), u(ni, nj), v(ni, nj), cp(ni), cp_exact(ni), q(ni), qa(ni))
  call circle_ring(a_cyl, ni, rx, ry)
  call build_ogrid(rx, ry, nj, router, 1.0e-3_real64, g, ierr)
  call grid_metrics(g, met)
  big_r = router
  k = big_r**2 / (big_r**2 - a_cyl**2)
  gamma_aux_exact = 2.0_real64 * pi / log(big_r / a_cyl)

  ! (a) prescribed c = 0: no circulation
  par%alpha_rad = 0.0_real64
  par%use_kutta = .false.
  par%c_prescribed = 0.0_real64
  par%verbose = .false.
  par%tol = 1.0e-13_real64
  call solve_flow(g, met, par, r0)
  call check(r0%ierr == ierr_ok, 'cylinder c=0 solve ok')
  print '(a,i0,a,i0,a,es9.2)', '  iterations stream ', r0%iters_stream, ' aux ', r0%iters_aux, ' resid ', r0%resid_stream
  call check_close(r0%gamma_stream_flux, 0.0_real64, tol_zero_cl, 'GammaStreamFlux = 0 for concentric cylinder')
  call velocity_field(met, r0%psi, u, v)
  call surface_cp(u, v, cp)
  itop = maxloc(ry, 1)
  top0 = hypot(u(itop, 1), v(itop, 1))
  do i = 1, ni
    th = atan2(ry(i) - center_y, rx(i) - center_x)
    speed_ex = 2.0_real64 * uinf * sin(th) * k
    cp_exact(i) = 1.0_real64 - speed_ex**2
  end do
  call check_close(maxval(abs(cp - cp_exact)), 0.0_real64, tol_cyl_cp_abs, 'surface Cp vs exact annulus (c = 0)')

  ! (b) measured GammaAux
  call check_close(r0%gamma_aux / gamma_aux_exact, 1.0_real64, tol_gammaaux_rel, 'GammaAux vs 2 pi / ln(R/a)')

  ! (c) prescribed circulation
  c = gamma_target * log(big_r / a_cyl) / (2.0_real64 * pi)
  par%c_prescribed = c
  call solve_flow(g, met, par, r1)
  call check(r1%ierr == ierr_ok, 'cylinder prescribed-circulation solve ok')
  call velocity_field(met, r1%psi, u, v)
  call surface_cp(u, v, cp)
  do i = 1, ni
    th = atan2(ry(i) - center_y, rx(i) - center_x)
    speed_ex = 2.0_real64 * uinf * sin(th) * k + gamma_target / (2.0_real64 * pi * a_cyl)
    cp_exact(i) = 1.0_real64 - speed_ex**2
  end do
  call check_close(maxval(abs(cp - cp_exact)), 0.0_real64, tol_cyl_cp_abs, 'surface Cp vs exact annulus (Gamma = 1)')
  gam = r1%gamma
  call check_close(2.0_real64 * gam / (uinf * 1.0_real64) / (2.0_real64 * gamma_target), 1.0_real64, tol_cyl_cl_rel, &
                   'CL = 2 Gamma / U vs exact')
  ! sign: clockwise Gamma raises speed at the top (largest y) node
  top1 = hypot(u(itop, 1), v(itop, 1))
  call check(top1 > top0, 'clockwise Gamma raises speed on upper side')
  call check(u(itop, 1) > 0.0_real64, 'flow at top of cylinder runs in +x')

  ! (d) Kutta denominator guard
  q = 1.0_real64
  qa = 0.0_real64
  par%use_kutta = .true.
  call kutta(g, q, qa, par, c, gam_aux, gam, ierr)
  call check(ierr == ierr_kutta, 'zero q_aux at Kutta nodes returns ierr_kutta')

  call finish()
end program test_cylinder
