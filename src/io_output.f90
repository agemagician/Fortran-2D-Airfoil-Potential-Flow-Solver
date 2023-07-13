! Run configuration (namelist), VTK legacy ASCII writer and surface Cp CSV writer.
module io_output
  use kinds, only: real64, ierr_ok, ierr_bad_input, ierr_io
  use grid, only: grid_record
  use naca_geometry, only: naca_ring
  use postproc, only: cp_field
  implicit none
  private
  public :: run_config, read_config, write_vtk, write_cp_csv

  type :: run_config
    character(len=16) :: naca_code = '2412'     ! longer than 4 characters is rejected
    real(real64) :: alpha_deg = 0.0_real64
    integer :: ni = 256, nj = 96
    real(real64) :: router = 30.0_real64
    real(real64) :: first_cell = 1.0e-3_real64
    real(real64) :: tol = 1.0e-10_real64
    integer :: maxiter = 2000
    character(len=256) :: vtk_path = 'flow.vtk'
    character(len=256) :: cp_path = 'cp.csv'
    character(len=256) :: input_path = 'input.nml'
  end type run_config

contains

  subroutine read_config(path, cfg, ierr, msg)
    character(*), intent(in) :: path
    type(run_config), intent(out) :: cfg
    integer, intent(out) :: ierr
    character(*), intent(out) :: msg
    character(len=16) :: naca_code
    real(real64) :: alpha_deg, router, first_cell, tol
    integer :: ni, nj, maxiter, u, ios, ie
    character(len=256) :: vtk_path, cp_path, iomsg
    real(real64) :: xt(32), yt(32)
    logical :: exists
    namelist /input/ naca_code, alpha_deg, ni, nj, router, first_cell, tol, maxiter, vtk_path, cp_path

    ierr = ierr_bad_input
    msg = ''
    naca_code = cfg%naca_code; alpha_deg = cfg%alpha_deg; ni = cfg%ni; nj = cfg%nj
    router = cfg%router; first_cell = cfg%first_cell; tol = cfg%tol; maxiter = cfg%maxiter
    vtk_path = cfg%vtk_path; cp_path = cfg%cp_path

    inquire (file=path, exist=exists)
    if (.not. exists) then
      msg = 'input file not found: '//trim(path)
      return
    end if
    open (newunit=u, file=path, status='old', action='read', iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
      msg = 'cannot open '//trim(path)//': '//trim(iomsg)
      return
    end if
    read (u, nml=input, iostat=ios, iomsg=iomsg)
    close (u)
    if (ios /= 0) then
      msg = 'namelist read error in '//trim(path)//': '//trim(iomsg)
      return
    end if

    if (len_trim(naca_code) /= 4) then
      msg = 'naca_code must be exactly four digits'
    else if (abs(alpha_deg) > 20.0_real64) then
      msg = 'alpha_deg must be between -20 and 20'
    else if (ni < 32 .or. mod(ni, 2) /= 0) then
      msg = 'ni must be even and at least 32'
    else if (nj < 16) then
      msg = 'nj must be at least 16'
    else if (router < 5.0_real64 .or. router > 100.0_real64) then
      msg = 'router must be between 5 and 100 chords'
    else if (trim(adjustl(vtk_path)) == trim(adjustl(cp_path))) then
      msg = 'vtk_path and cp_path must differ'
    else if (trim(adjustl(vtk_path)) == trim(adjustl(path)) .or. trim(adjustl(cp_path)) == trim(adjustl(path))) then
      msg = 'an output path equals the input path'
    end if
    if (len_trim(msg) == 0) then
      call naca_ring(trim(naca_code), 32, xt, yt, ie)
      if (ie /= ierr_ok) msg = 'naca_code is not a valid four-digit code with nonzero thickness'
    end if
    if (len_trim(msg) > 0) return

    cfg%naca_code = naca_code; cfg%alpha_deg = alpha_deg; cfg%ni = ni; cfg%nj = nj
    cfg%router = router; cfg%first_cell = first_cell; cfg%tol = tol; cfg%maxiter = maxiter
    cfg%vtk_path = vtk_path; cfg%cp_path = cp_path; cfg%input_path = path
    ierr = ierr_ok
  end subroutine read_config

  subroutine write_vtk(path, input_path, g, psi, u, v, ierr)
    character(*), intent(in) :: path, input_path
    type(grid_record), intent(in) :: g
    real(real64), intent(in) :: psi(:, :), u(:, :), v(:, :)
    integer, intent(out) :: ierr
    integer :: unit, ios, i, j, ii, npts
    character(len=256) :: iomsg

    ierr = ierr_ok
    if (trim(adjustl(path)) == trim(adjustl(input_path))) then
      print '(a)', 'error: refusing to overwrite the input file with '//trim(path)
      ierr = ierr_io
      return
    end if
    open (newunit=unit, file=path, status='replace', action='write', iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
      print '(a)', 'error: cannot write '//trim(path)//': '//trim(iomsg)
      ierr = ierr_io
      return
    end if
    npts = (g%ni + 1) * g%nj
    write (unit, '(a)') '# vtk DataFile Version 3.0'
    write (unit, '(a)') 'airfoil potential flow'
    write (unit, '(a)') 'ASCII'
    write (unit, '(a)') 'DATASET STRUCTURED_GRID'
    write (unit, '(a,i0,1x,i0,a)') 'DIMENSIONS ', g%ni + 1, g%nj, ' 1'
    write (unit, '(a,i0,a)') 'POINTS ', npts, ' double'
    do j = 1, g%nj
      do i = 1, g%ni + 1
        ii = mod(i - 1, g%ni) + 1
        write (unit, '(3es24.15e3)') g%x(ii, j), g%y(ii, j), 0.0_real64
      end do
    end do
    write (unit, '(a,i0)') 'POINT_DATA ', npts
    write (unit, '(a)') 'SCALARS psi double 1'
    write (unit, '(a)') 'LOOKUP_TABLE default'
    do j = 1, g%nj
      do i = 1, g%ni + 1
        write (unit, '(es24.15e3)') psi(mod(i - 1, g%ni) + 1, j)
      end do
    end do
    write (unit, '(a)') 'SCALARS cp double 1'
    write (unit, '(a)') 'LOOKUP_TABLE default'
    do j = 1, g%nj
      do i = 1, g%ni + 1
        ii = mod(i - 1, g%ni) + 1
        write (unit, '(es24.15e3)') cp_field(u(ii, j), v(ii, j))
      end do
    end do
    write (unit, '(a)') 'VECTORS velocity double'
    do j = 1, g%nj
      do i = 1, g%ni + 1
        ii = mod(i - 1, g%ni) + 1
        write (unit, '(3es24.15e3)') u(ii, j), v(ii, j), 0.0_real64
      end do
    end do
    close (unit, iostat=ios)
    if (ios /= 0) then
      print '(a)', 'error: write failed for '//trim(path)
      ierr = ierr_io
    end if
  end subroutine write_vtk

  subroutine write_cp_csv(path, input_path, g, cp_ring, ierr)
    character(*), intent(in) :: path, input_path
    type(grid_record), intent(in) :: g
    real(real64), intent(in) :: cp_ring(:)
    integer, intent(out) :: ierr
    integer :: unit, ios, i
    character(len=256) :: iomsg

    ierr = ierr_ok
    if (trim(adjustl(path)) == trim(adjustl(input_path))) then
      print '(a)', 'error: refusing to overwrite the input file with '//trim(path)
      ierr = ierr_io
      return
    end if
    open (newunit=unit, file=path, status='replace', action='write', iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
      print '(a)', 'error: cannot write '//trim(path)//': '//trim(iomsg)
      ierr = ierr_io
      return
    end if
    write (unit, '(a)') 'x,y,cp'
    do i = 1, g%ni
      write (unit, '(es23.15e3,a,es23.15e3,a,es23.15e3)') g%x(i, 1), ',', g%y(i, 1), ',', cp_ring(i)
    end do
    close (unit, iostat=ios)
    if (ios /= 0) then
      print '(a)', 'error: write failed for '//trim(path)
      ierr = ierr_io
    end if
  end subroutine write_cp_csv

end module io_output
