program test_grid_diagnostics
  use kinds
  use constants
  use naca_geometry
  use grid
  use reference_geometry
  use csr_solver
  use potential_flow
  use tolerances
  use test_support
  implicit none

  type(grid_record) :: g, g2
  type(quality_report) :: rep
  real(real64), allocatable :: x(:), y(:)
  real(real64) :: asum, expected
  type(metrics_record) :: met
  type(csr_matrix) :: a
  type(solver_params) :: par
  type(flow_result) :: r_a, r_b
  real(real64), allocatable :: psi_e(:, :), vec(:), res(:), bvec(:)
  logical :: sorted
  integer :: ierr, i, j, ni, nj, k, s_shift

  ! default case: NACA 2412, 256 x 96, 30 chords, first cell 1e-3
  ni = 256; nj = 96
  allocate (x(ni), y(ni))
  call naca_ring('2412', ni, x, y, ierr)
  call build_ogrid(x, y, nj, 30.0_real64, 1.0e-3_real64, g, ierr)
  call check(ierr == ierr_ok, 'default grid builds (spacing feasible)')

  asum = 0.0_real64
  do j = 1, nj - 1
    do i = 1, ni
      asum = asum + cell_area(g, i, j)
    end do
  end do
  expected = abs(shoelace(g%x(:, nj), g%y(:, nj))) - abs(shoelace(g%x(:, 1), g%y(:, 1)))
  call check_close(asum / expected, 1.0_real64, 1.0e-12_real64, 'cell areas sum to outer minus body polygon')
  call check(cell_area(g, ni, 1) > 0.0_real64, 'seam cell (nodes ni and 1) has positive area')

  call check_quality(g, rep, ierr)
  print '(a,es10.3,a,f7.2,a,es10.3,a,es10.3,a,es10.3)', '  min_area ', rep%min_area, '  skew ', rep%max_skew_deg, &
    '  aspect ', rep%max_aspect, '  area_ratio ', rep%area_ratio, '  min_jac ', rep%min_jac
  call check(ierr == ierr_ok, 'default grid passes the quality gate')
  call check(rep%max_aspect * 2.0_real64 <= max_aspect, 'aspect margin factor 2')
  call check(rep%area_ratio >= 2.0_real64 * min_area_ratio, 'area ratio margin factor 2')
  call check(rep%max_skew_deg <= max_skew_deg - 5.0_real64, 'skew margin 5 degrees')

  ! folded grid: interior node moved past its radial neighbour
  g2 = g
  g2%x(100, 50) = g%x(100, 52)
  g2%y(100, 50) = g%y(100, 52)
  call check_quality(g2, rep, ierr)
  call check(ierr == ierr_quality .and. trim(rep%limit) == 'min_area', 'folded grid rejected, names min_area')
  call check(rep%first_bad_j >= 48 .and. rep%first_bad_j <= 51, 'first offending cell located')

  ! limits forced past: report values -> limit name
  rep = quality_report(1.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64, 0, 0, '')
  call check(len_trim(quality_failure(rep)) == 0, 'clean report passes')
  rep%max_skew_deg = max_skew_deg + 1.0_real64
  call check(trim(quality_failure(rep)) == 'max_skew', 'skew limit named')
  rep%max_skew_deg = 0.0_real64; rep%max_aspect = max_aspect * 2.0_real64
  call check(trim(quality_failure(rep)) == 'max_aspect', 'aspect limit named')
  rep%max_aspect = 1.0_real64; rep%area_ratio = min_area_ratio / 2.0_real64
  call check(trim(quality_failure(rep)) == 'area_ratio', 'area ratio limit named')
  rep%area_ratio = 1.0_real64; rep%min_jac = -1.0_real64
  call check(trim(quality_failure(rep)) == 'min_jacobian', 'nodal Jacobian limit named')

  ! degenerate cell: node pushed onto its tangential neighbour
  g2 = g
  g2%x(100, 40) = g%x(101, 40)
  g2%y(100, 40) = g%y(101, 40)
  call check_quality(g2, rep, ierr)
  call check(ierr == ierr_quality, 'collapsed edge rejected by the gate')
  print '(a,a)', '  sheared grid limit: ', trim(rep%limit)

  ! tiny first cell: aspect ratio limit
  call build_ogrid(x, y, 24, 30.0_real64, 1.0e-6_real64, g2, ierr)
  call check_quality(g2, rep, ierr)
  print '(a,a,es10.3,es10.3)', '  tiny first_cell limit: ', trim(rep%limit), rep%max_aspect, rep%area_ratio
  call check(ierr == ierr_quality .and. trim(rep%limit) == 'max_aspect', 'tiny first_cell rejected on aspect ratio')

  ! infeasible first_cell
  call build_ogrid(x, y, nj, 30.0_real64, 0.5_real64, g2, ierr)
  call check(ierr == ierr_spacing, 'infeasible first_cell returns ierr_spacing')

  ! ---- phase 5 additions ----
  ! CSR pattern: sorted columns, nnz formula
  call grid_metrics(g, met)
  call assemble(g, met, a)
  sorted = .true.
  do i = 1, a%n
    do k = a%rowptr(i) + 1, a%rowptr(i + 1) - 1
      if (a%colidx(k) <= a%colidx(k - 1)) sorted = .false.
    end do
  end do
  call check(sorted, 'CSR columns sorted in every row (wrap rows included)')
  call check(a%nnz == 9 * ni * (nj - 2) - 6 * ni .and. a%rowptr(a%n + 1) - 1 == a%nnz, 'CSR nnz matches formula')

  ! linear-psi patch test on a uniform parallelogram grid (rows away from the periodic seam)
  deallocate (g2%x, g2%y)
  g2%ni = 64; g2%nj = 20
  allocate (g2%x(64, 20), g2%y(64, 20))
  do j = 1, 20
    do i = 1, 64
      g2%x(i, j) = 0.10_real64 * i + 0.03_real64 * j
      g2%y(i, j) = 0.02_real64 * i + 0.08_real64 * j
    end do
  end do
  call grid_metrics(g2, met)
  call assemble(g2, met, a)
  allocate (psi_e(64, 20), vec(a%n), res(a%n), bvec(a%n))
  psi_e = 1.3_real64 * g2%x - 0.7_real64 * g2%y
  do j = 2, 19
    vec((j - 2) * 64 + 1:(j - 1) * 64) = psi_e(:, j)
  end do
  call csr_matvec(a, vec, res)
  call build_rhs_dirichlet(g2, met, psi_e(:, 1), psi_e(:, 20), bvec)
  res = res - bvec
  call check_close(patch_max(res), 0.0_real64, tol_patch, 'linear psi patch test (stencil + ring terms)')

  ! rotated-seam test: same grid, i origin shifted; CL must not change
  par%alpha_rad = 5.0_real64 * pi / 180.0_real64
  par%tol = 1.0e-12_real64
  par%maxiter = 5000
  call naca_ring('2412', 128, x(1:128), y(1:128), ierr)
  call build_ogrid(x(1:128), y(1:128), 48, 30.0_real64, 2.0e-3_real64, g, ierr)
  call grid_metrics(g, met)
  call solve_flow(g, met, par, r_a)
  s_shift = 37
  g2%ni = g%ni; g2%nj = g%nj
  deallocate (g2%x, g2%y)
  allocate (g2%x(g%ni, g%nj), g2%y(g%ni, g%nj))
  do i = 1, g%ni
    g2%x(i, :) = g%x(mod(i - 1 + s_shift, g%ni) + 1, :)
    g2%y(i, :) = g%y(mod(i - 1 + s_shift, g%ni) + 1, :)
  end do
  call grid_metrics(g2, met)
  par%kutta_i_lower = modulo(2 - 1 - s_shift, g%ni) + 1
  par%kutta_i_upper = modulo(g%ni - 1 - 1 - s_shift, g%ni) + 1
  call solve_flow(g2, met, par, r_b)
  print '(a,es14.6,a,es14.6)', '  CL unshifted ', 2.0_real64 * r_a%gamma, '  shifted ', 2.0_real64 * r_b%gamma
  call check(r_a%ierr == ierr_ok .and. r_b%ierr == ierr_ok, 'rotated-seam solves converge')
  call check_close(2.0_real64 * r_b%gamma, 2.0_real64 * r_a%gamma, tol_rotated_seam_cl, 'CL invariant under seam rotation')

  call finish()

contains

  ! largest residual among rows whose stencil does not touch the periodic seam (i = 3 .. 62)
  function patch_max(r) result(m)
    real(real64), intent(in) :: r(:)
    real(real64) :: m
    integer :: jj, ii
    m = 0.0_real64
    do jj = 2, 19
      do ii = 3, 62
        m = max(m, abs(r((jj - 2) * 64 + ii)))
      end do
    end do
  end function patch_max
end program test_grid_diagnostics
