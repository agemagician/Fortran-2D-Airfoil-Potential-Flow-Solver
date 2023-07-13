! Run by `make convergence`, not by `make test`. Prints the study that sets tests/tolerances.f90.
program test_grid_convergence
  use kinds
  use constants, only: pi
  use grid
  use potential_flow
  use postproc
  use naca_geometry
  use reference_geometry
  use tolerances
  use test_support
  implicit none

  integer, parameter :: levels(2, 3) = reshape([128, 48, 256, 96, 512, 192], [2, 3])
  real(real64), parameter :: a_cyl = 0.5_real64, gamma_target = 1.0_real64
  real(real64) :: cyl_err(3), cyl_gam(3), jk_err(3), jk_gam(3), cl, cl_exact, change_cyl, change_jk
  real(real64) :: cl_r(3), routers(3), e(3), cl_e(3)
  integer :: l, k

  print '(a)', 'cylinder (prescribed Gamma = 1)'
  do l = 1, 3
    call cylinder_case(levels(1, l), levels(2, l), cl, cyl_gam(l))
    cyl_err(l) = abs(cl - 2.0_real64 * gamma_target)
    call report(l, cyl_err, cyl_gam)
  end do
  change_cyl = abs(cyl_gam(3) - cyl_gam(2)) / abs(cyl_gam(3))

  print '(a)', 'Joukowski blunt-face limit (e = 6.25e-4, alpha = 4 deg)'
  cl_exact = joukowski_cl(4.0_real64 * pi / 180.0_real64)
  do l = 1, 3
    call joukowski_case(6.25e-4_real64, levels(1, l), levels(2, l), 30.0_real64, cl, jk_gam(l))
    jk_err(l) = abs(cl - cl_exact)
    call report(l, jk_err, jk_gam)
  end do
  change_jk = abs(jk_gam(3) - jk_gam(2)) / abs(jk_gam(3))

  print '(a)', 'router sweep, NACA 2412 at 5 deg, 256 x 96'
  routers = [20.0_real64, 30.0_real64, 50.0_real64]
  do k = 1, 3
    call naca_case(routers(k), cl_r(k))
  end do
  do k = 1, 3
    print '(a,f6.1,a,f10.6,a,es10.3)', '  router ', routers(k), '  CL ', cl_r(k), '  shift vs 30: ', cl_r(k) - cl_r(2)
  end do

  print '(a)', 'lift against blunt-face thickness, Joukowski at 256 x 96'
  e = [2.5e-3_real64, 1.25e-3_real64, 6.25e-4_real64]
  do k = 1, 3
    call joukowski_case(e(k), 256, 96, 30.0_real64, cl_e(k), cl)
    print '(a,es10.3,a,f10.6,a,es10.3)', '  e ', e(k), '  CL ', cl_e(k), '  error ', cl_e(k) - cl_exact
  end do

  print '(a,es10.3,a,es10.3)', 'Gamma change, two finest levels: cylinder ', change_cyl, '  Joukowski ', change_jk
  call check(change_cyl <= tol_gamma_change_rel, 'cylinder Gamma change within tolerance')
  call check(change_jk <= tol_gamma_change_rel, 'Joukowski Gamma change within tolerance')
  call finish()

contains

  subroutine report(l_, err, gam)
    integer, intent(in) :: l_
    real(real64), intent(in) :: err(3), gam(3)
    real(real64) :: chg, ord
    chg = 0.0_real64
    ord = 0.0_real64
    if (l_ > 1) then
      chg = abs(gam(l_) - gam(l_ - 1)) / abs(gam(l_))
      if (err(l_) > 0.0_real64) ord = log(err(l_ - 1) / err(l_)) / log(2.0_real64)
    end if
    print '(a,i5,a,i4,a,es10.3,a,es10.3,a,f6.2)', '  ni ', levels(1, l_), ' nj ', levels(2, l_), '  CL error ', err(l_), &
      '  Gamma change ', chg, '  order ', ord
  end subroutine report

  subroutine solve(rx, ry, njl, router, alpha_rad, par, r)
    real(real64), intent(in) :: rx(:), ry(:), router, alpha_rad
    integer, intent(in) :: njl
    type(solver_params), intent(inout) :: par
    type(flow_result), intent(out) :: r
    type(grid_record) :: g
    type(metrics_record) :: met
    integer :: ierr_l
    call build_ogrid(rx, ry, njl, router, 1.0e-3_real64, g, ierr_l)
    call grid_metrics(g, met)
    par%alpha_rad = alpha_rad
    par%tol = 1.0e-10_real64
    par%maxiter = 5000
    call solve_flow(g, met, par, r)
  end subroutine solve

  subroutine cylinder_case(nic, njc, cl_out, gam_out)
    integer, intent(in) :: nic, njc
    real(real64), intent(out) :: cl_out, gam_out
    real(real64) :: rx(nic), ry(nic)
    type(solver_params) :: par
    type(flow_result) :: r
    call circle_ring(a_cyl, nic, rx, ry)
    par%use_kutta = .false.
    par%c_prescribed = gamma_target * log(30.0_real64 / a_cyl) / (2.0_real64 * pi)
    call solve(rx, ry, njc, 30.0_real64, 0.0_real64, par, r)
    cl_out = 2.0_real64 * r%gamma
    gam_out = r%gamma
  end subroutine cylinder_case

  subroutine joukowski_case(et, nic, njc, router, cl_out, gam_out)
    real(real64), intent(in) :: et, router
    integer, intent(in) :: nic, njc
    real(real64), intent(out) :: cl_out, gam_out
    real(real64) :: rx(nic), ry(nic)
    type(solver_params) :: par
    type(flow_result) :: r
    call joukowski_ring(et, nic, rx, ry)
    call solve(rx, ry, njc, router, 4.0_real64 * pi / 180.0_real64, par, r)
    cl_out = 2.0_real64 * r%gamma
    gam_out = r%gamma
  end subroutine joukowski_case

  subroutine naca_case(router, cl_out)
    real(real64), intent(in) :: router
    real(real64), intent(out) :: cl_out
    real(real64) :: rx(256), ry(256)
    type(solver_params) :: par
    type(flow_result) :: r
    integer :: ierr_l
    call naca_ring('2412', 256, rx, ry, ierr_l)
    call solve(rx, ry, 96, router, 5.0_real64 * pi / 180.0_real64, par, r)
    cl_out = 2.0_real64 * r%gamma
  end subroutine naca_case

end program test_grid_convergence
