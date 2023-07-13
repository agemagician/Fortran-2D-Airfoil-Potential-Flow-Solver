program airfoil
  use kinds
  use constants, only: pi, cl_gap_flag
  use naca_geometry, only: naca_ring
  use grid
  use potential_flow
  use postproc
  use io_output
  implicit none

  type(run_config) :: cfg
  type(grid_record) :: g
  type(metrics_record) :: met
  type(quality_report) :: rep
  type(solver_params) :: par
  type(flow_result) :: res
  type(force_report) :: fr
  real(real64), allocatable :: x(:), y(:), u(:, :), v(:, :), cp(:)
  character(len=256) :: path, msg
  integer :: ierr, nargs
  real(real64) :: cl_circ

  path = 'input.nml'
  nargs = command_argument_count()
  if (nargs >= 1) call get_command_argument(1, path)

  call read_config(trim(path), cfg, ierr, msg)
  if (ierr /= ierr_ok) call fail('config', trim(msg))

  allocate (x(cfg%ni), y(cfg%ni))
  call naca_ring(trim(cfg%naca_code), cfg%ni, x, y, ierr)
  if (ierr /= ierr_ok) call fail('geometry', 'cannot build the airfoil ring')

  call build_ogrid(x, y, cfg%nj, cfg%router, cfg%first_cell, g, ierr)
  if (ierr /= ierr_ok) call fail('grid', 'cannot build the O-grid (see message above)')
  call check_quality(g, rep, ierr)
  print '(a,es12.5)', 'grid min_area: ', rep%min_area
  print '(a,es12.5)', 'grid max_skew_deg: ', rep%max_skew_deg
  print '(a,es12.5)', 'grid max_aspect: ', rep%max_aspect
  print '(a,es12.5)', 'grid area_ratio: ', rep%area_ratio
  if (ierr /= ierr_ok) then
    write (msg, '(a,a,a,i0,a,i0,a)') 'limit ', trim(rep%limit), ' failed, first offending cell (i=', &
      rep%first_bad_i, ', j=', rep%first_bad_j, ')'
    call fail('grid quality', trim(msg))
  end if
  call grid_metrics(g, met)

  par%alpha_rad = cfg%alpha_deg * pi / 180.0_real64
  par%tol = cfg%tol
  par%maxiter = cfg%maxiter
  par%verbose = .true.
  call solve_flow(g, met, par, res)
  select case (res%ierr)
  case (ierr_ok)
  case (ierr_pivot)
    call fail('solve', 'ILU(0) pivot too small')
  case (ierr_breakdown)
    call fail('solve', 'BiCGSTAB breakdown')
  case (ierr_noconv)
    write (msg, '(a,es10.3,a,es10.3)') 'BiCGSTAB did not converge; residuals ', res%resid_stream, ' / ', res%resid_aux
    call fail('solve', trim(msg))
  case (ierr_kutta)
    call fail('solve', 'Kutta condition denominator too small')
  case default
    call fail('solve', 'solver failed')
  end select

  allocate (u(cfg%ni, cfg%nj), v(cfg%ni, cfg%nj), cp(cfg%ni))
  call velocity_field(met, res%psi, u, v)
  call surface_cp(u, v, cp)
  cl_circ = 2.0_real64 * res%gamma
  call pressure_forces(g, cp, par%alpha_rad, cl_circ, fr)

  print '(a,i0)', 'iters_stream: ', res%iters_stream
  print '(a,es12.5)', 'resid_stream: ', res%resid_stream
  print '(a,i0)', 'iters_aux: ', res%iters_aux
  print '(a,es12.5)', 'resid_aux: ', res%resid_aux
  print '(a,es12.5)', 'gamma_stream_flux: ', res%gamma_stream_flux
  print '(a,es12.5)', 'gamma: ', res%gamma
  print '(a,es12.5)', 'cl_circ: ', cl_circ
  print '(a,es12.5)', 'cl_press: ', fr%cl_press
  print '(a,es12.5)', 'cl_gap: ', fr%cl_gap_rel
  print '(a,es12.5)', 'cd: ', fr%cd
  if (fr%cl_gap_rel > cl_gap_flag) print '(a,es12.5,a,es12.5)', 'FLAG: cl gap ', fr%cl_gap_rel, ' exceeds ', cl_gap_flag

  call write_vtk(trim(cfg%vtk_path), trim(cfg%input_path), g, res%psi, u, v, ierr)
  if (ierr /= ierr_ok) error stop 1
  call write_cp_csv(trim(cfg%cp_path), trim(cfg%input_path), g, cp, ierr)
  if (ierr /= ierr_ok) error stop 1

contains

  subroutine fail(step, text)
    character(*), intent(in) :: step, text
    print '(a)', 'error ['//step//']: '//text
    error stop 1
  end subroutine fail

end program airfoil
