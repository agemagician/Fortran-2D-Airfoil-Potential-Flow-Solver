program test_io_and_config
  use kinds
  use grid
  use io_output
  use reference_geometry
  use test_support
  implicit none

  character(len=*), parameter :: dir = 'build/debug/'
  type(run_config) :: cfg
  type(grid_record) :: g
  real(real64), allocatable :: rx(:), ry(:), psi(:, :), u(:, :), v(:, :), cp(:)
  character(len=256) :: msg, line, first_pt, seam_pt
  integer :: ierr, ni, nj, unit, ios, i, nrows, ierr2

  ! --- read_config ---
  call put_nml('t_good.nml', "naca_code='0012', alpha_deg=3.0, ni=64, nj=24, vtk_path='build/debug/o.vtk', cp_path='build/debug/o.csv'")
  call read_config(dir//'t_good.nml', cfg, ierr, msg)
  call check(ierr == ierr_ok .and. cfg%ni == 64 .and. trim(cfg%naca_code) == '0012', 'valid namelist read')

  call read_config(dir//'missing.nml', cfg, ierr, msg)
  call check(ierr == ierr_bad_input, 'missing file rejected')
  open (newunit=unit, file=dir//'t_bad.nml', status='replace')
  write (unit, '(a)') 'this is not a namelist'
  close (unit)
  call read_config(dir//'t_bad.nml', cfg, ierr, msg)
  call check(ierr == ierr_bad_input, 'malformed namelist rejected')

  call expect_bad("naca_code='24'", 'short naca_code')
  call expect_bad("naca_code='24123'", 'long naca_code')
  call expect_bad("naca_code='2400'", 'zero thickness code')
  call expect_bad("alpha_deg=25.0", 'alpha out of range')
  call expect_bad("ni=7", 'ni too small')
  call expect_bad("ni=33", 'ni odd')
  call expect_bad("nj=8", 'nj too small')
  call expect_bad("router=200.0", 'router out of range')
  call expect_bad("vtk_path='a.out', cp_path='a.out'", 'vtk_path equals cp_path')
  call expect_bad("vtk_path='build/debug/t_x.nml'", 'output path equals input path')

  ! --- writers on a small circle grid ---
  ni = 32; nj = 16
  allocate (rx(ni), ry(ni), psi(ni, nj), u(ni, nj), v(ni, nj), cp(ni))
  call circle_ring(0.5_real64, ni, rx, ry)
  call build_ogrid(rx, ry, nj, 10.0_real64, 1.0e-2_real64, g, ierr)
  psi = g%y
  u = 1.0_real64
  v = 0.0_real64
  cp = 0.25_real64
  call write_vtk(dir//'t.vtk', dir//'t_in.nml', g, psi, u, v, ierr)
  call check(ierr == ierr_ok, 'write_vtk ok')
  open (newunit=unit, file=dir//'t.vtk', status='old', action='read')
  read (unit, '(a)') line
  call check(trim(line) == '# vtk DataFile Version 3.0', 'vtk first header line')
  do i = 1, 3
    read (unit, '(a)') line
  end do
  call check(trim(line) == 'DATASET STRUCTURED_GRID', 'vtk dataset line')
  read (unit, '(a)') line
  call check(trim(line) == 'DIMENSIONS 33 16 1', 'vtk DIMENSIONS line')
  read (unit, '(a)') line
  call check(trim(line) == 'POINTS 528 double', 'vtk point count')
  read (unit, '(a)') first_pt
  do i = 2, 33
    read (unit, '(a)') seam_pt
  end do
  call check(trim(first_pt) == trim(seam_pt), 'seam column repeated at i = ni + 1')
  nrows = 0
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:10) == 'POINT_DATA' .or. line(1:7) == 'SCALARS' .or. line(1:7) == 'VECTORS') nrows = nrows + 1
  end do
  close (unit)
  call check(nrows == 4, 'POINT_DATA, two SCALARS and VECTORS present')

  call write_cp_csv(dir//'t.csv', dir//'t_in.nml', g, cp, ierr)
  call check(ierr == ierr_ok, 'write_cp_csv ok')
  open (newunit=unit, file=dir//'t.csv', status='old', action='read')
  read (unit, '(a)') line
  call check(trim(line) == 'x,y,cp', 'csv header')
  nrows = 0
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    nrows = nrows + 1
  end do
  close (unit)
  call check(nrows == ni, 'csv has ni rows')

  call write_vtk(dir//'t_in.nml', dir//'t_in.nml', g, psi, u, v, ierr)
  call write_cp_csv(dir//'t_in.nml', dir//'t_in.nml', g, cp, ierr2)
  call check(ierr == ierr_io .and. ierr2 == ierr_io, 'writers refuse to overwrite the input path')

  call finish()

contains

  subroutine put_nml(name, body)
    character(*), intent(in) :: name, body
    integer :: un
    open (newunit=un, file=dir//name, status='replace')
    write (un, '(a)') '&input'
    write (un, '(a)') body
    write (un, '(a)') '/'
    close (un)
  end subroutine put_nml

  subroutine expect_bad(body, label)
    character(*), intent(in) :: body, label
    type(run_config) :: c
    character(len=256) :: m
    integer :: ie
    call put_nml('t_x.nml', body)
    call read_config(dir//'t_x.nml', c, ie, m)
    call check(ie == ierr_bad_input, label//' rejected')
  end subroutine expect_bad

end program test_io_and_config
