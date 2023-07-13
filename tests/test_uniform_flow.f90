program test_uniform_flow
  use kinds
  use constants, only: uinf
  use grid
  use csr_solver
  use potential_flow
  use postproc
  use reference_geometry
  use tolerances
  use test_support
  implicit none

  real(real64), parameter :: alpha = 0.3_real64
  type(grid_record) :: g
  type(metrics_record) :: met
  real(real64), allocatable :: u(:, :), v(:, :), psi(:, :)
  real(real64) :: e_coarse(2), e_fine(2)
  integer :: i, j, ni, nj

  ! (a) uniform parallelogram grid: exact psi gives exact (u, v) away from the periodic seam
  ni = 40; nj = 20
  g%ni = ni; g%nj = nj
  allocate (g%x(ni, nj), g%y(ni, nj), psi(ni, nj), u(ni, nj), v(ni, nj))
  do j = 1, nj
    do i = 1, ni
      g%x(i, j) = 0.10_real64 * i + 0.03_real64 * j
      g%y(i, j) = 0.02_real64 * i + 0.08_real64 * j
    end do
  end do
  psi = uinf * (g%y * cos(alpha) - g%x * sin(alpha))
  call grid_metrics(g, met)
  call velocity_field(met, psi, u, v)
  call check_close(maxval(abs(u(2:ni - 1, :) - uinf * cos(alpha))), 0.0_real64, tol_velocity_exact, 'u = U cos(alpha)')
  call check_close(maxval(abs(v(2:ni - 1, :) - uinf * sin(alpha))), 0.0_real64, tol_velocity_exact, 'v = U sin(alpha)')

  ! (b) circular grid, uniform stream as Dirichlet data on both rings
  call solve_uniform(128, 48, 2.0e-3_real64, e_coarse)
  call solve_uniform(256, 96, 1.0e-3_real64, e_fine)
  print '(a,2es11.3,a,2es11.3)', '  psi,vel error 128x48: ', e_coarse, '   256x96: ', e_fine
  call check(e_coarse(1) / e_fine(1) >= uniform_refine_ratio_min, 'psi error falls with refinement')
  call check(e_coarse(2) / e_fine(2) >= uniform_refine_ratio_min, 'velocity error falls with refinement')

  call finish()

contains

  subroutine solve_uniform(nic, njc, h, errs)
    integer, intent(in) :: nic, njc
    real(real64), intent(in) :: h
    real(real64), intent(out) :: errs(2)
    type(grid_record) :: gc
    type(metrics_record) :: mc
    type(csr_matrix) :: a, lu
    real(real64), allocatable :: rx(:), ry(:), b(:), x(:), p(:, :), pe(:, :), uu(:, :), vv(:, :)
    real(real64) :: resid
    integer :: ierr, it, jj

    allocate (rx(nic), ry(nic))
    call circle_ring(0.5_real64, nic, rx, ry)
    call build_ogrid(rx, ry, njc, 30.0_real64, h, gc, ierr)
    call grid_metrics(gc, mc)
    call assemble(gc, mc, a)
    call ilu0_factor(a, lu, ierr)
    allocate (b(a%n), x(a%n), p(nic, njc), pe(nic, njc), uu(nic, njc), vv(nic, njc))
    pe = uinf * (gc%y * cos(alpha) - gc%x * sin(alpha))
    call build_rhs_dirichlet(gc, mc, pe(:, 1), pe(:, njc), b)
    call bicgstab(a, lu, b, 1.0e-12_real64, 3000, .false., x, it, resid, ierr)
    call check(ierr == ierr_ok, 'uniform-flow solve converges')
    print '(a,i0,a,i0,a,i0,a,es9.2)', '  ', nic, 'x', njc, ' iterations ', it, ' resid ', resid
    p(:, 1) = pe(:, 1)
    p(:, njc) = pe(:, njc)
    do jj = 2, njc - 1
      p(:, jj) = x((jj - 2) * nic + 1:(jj - 1) * nic)
    end do
    call velocity_field(mc, p, uu, vv)
    errs(1) = maxval(abs(p - pe))
    errs(2) = maxval(abs(uu(:, 1:njc - 1) - uinf * cos(alpha))) + maxval(abs(vv(:, 1:njc - 1) - uinf * sin(alpha)))
  end subroutine solve_uniform

end program test_uniform_flow
