program test_joukowski
  use kinds
  use constants, only: pi
  use naca_geometry
  use grid
  use potential_flow
  use reference_geometry
  use tolerances
  use test_support
  implicit none

  integer, parameter :: ni = 256, nj = 96
  real(real64), parameter :: alpha = 4.0_real64 * pi / 180.0_real64
  real(real64) :: x(ni), y(ni), cl(3), e(3), cl_exact, err1, err2, cinf
  integer :: k, ierr

  cl_exact = joukowski_cl(alpha)
  print '(a,f9.5)', '  exact sharp-edge CL at 4 deg: ', cl_exact
  e = [2.5e-3_real64, 1.25e-3_real64, 6.25e-4_real64]
  do k = 1, 3
    call joukowski_ring(e(k), ni, x, y)
    call circulation_cl(x, y, alpha, cl(k), ierr)
    call check(ierr == ierr_ok, 'joukowski solve ok')
    print '(a,es10.3,a,f9.5,a,es10.3)', '  e = ', e(k), '  CL = ', cl(k), '  error = ', cl(k) - cl_exact
  end do
  err1 = abs(cl(1) - cl_exact)
  err2 = abs(cl(2) - cl_exact)
  call check(err2 <= joukowski_ratio_max * err1, 'CL error falls when the face thickness halves')
  cinf = cl(3) - (cl(3) - cl(2))**2 / ((cl(3) - cl(2)) - (cl(2) - cl(1)))
  print '(a,f9.5)', '  Aitken extrapolated CL: ', cinf
  call check_close(cinf / cl_exact, 1.0_real64, tol_joukowski_extrap_rel, 'Aitken limit matches exact Joukowski CL')

  ! symmetric NACA 0012 at alpha = 0: Kutta path on a blunt edge with an exact answer
  call naca_ring('0012', ni, x, y, ierr)
  call circulation_cl(x, y, 0.0_real64, cl(1), ierr)
  call check_close(cl(1), 0.0_real64, tol_zero_cl, 'NACA 0012 alpha=0 CL from circulation')

  call finish()

contains

  subroutine circulation_cl(rx, ry, alpha_rad, cl_out, ierr_out)
    real(real64), intent(in) :: rx(:), ry(:), alpha_rad
    real(real64), intent(out) :: cl_out
    integer, intent(out) :: ierr_out
    type(grid_record) :: g
    type(metrics_record) :: met
    type(solver_params) :: par
    type(flow_result) :: r
    call build_ogrid(rx, ry, nj, 30.0_real64, 1.0e-3_real64, g, ierr_out)
    if (ierr_out /= ierr_ok) return
    call grid_metrics(g, met)
    par%alpha_rad = alpha_rad
    par%maxiter = 5000
    par%tol = 1.0e-13_real64
    call solve_flow(g, met, par, r)
    ierr_out = r%ierr
    cl_out = 2.0_real64 * r%gamma
  end subroutine circulation_cl

end program test_joukowski
