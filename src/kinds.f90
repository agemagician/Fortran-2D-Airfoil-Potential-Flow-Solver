module kinds
  implicit none
  private
  public :: real64
  public :: ierr_ok, ierr_bad_input, ierr_spacing, ierr_quality, ierr_pivot
  public :: ierr_breakdown, ierr_noconv, ierr_kutta, ierr_io

  integer, parameter :: real64 = selected_real_kind(15, 307)

  integer, parameter :: ierr_ok = 0
  integer, parameter :: ierr_bad_input = 1
  integer, parameter :: ierr_spacing = 2
  integer, parameter :: ierr_quality = 3
  integer, parameter :: ierr_pivot = 4
  integer, parameter :: ierr_breakdown = 5
  integer, parameter :: ierr_noconv = 6
  integer, parameter :: ierr_kutta = 7
  integer, parameter :: ierr_io = 8
end module kinds
