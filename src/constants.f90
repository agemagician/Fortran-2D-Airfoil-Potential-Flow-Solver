module constants
  use kinds, only: real64
  implicit none
  public

  real(real64), parameter :: pi = 3.14159265358979323846_real64
  real(real64), parameter :: uinf = 1.0_real64
  real(real64), parameter :: chord = 1.0_real64
  real(real64), parameter :: center_x = 0.25_real64
  real(real64), parameter :: center_y = 0.0_real64

  ! grid quality limits (min area and min Jacobian must exceed zero)
  real(real64), parameter :: max_skew_deg = 80.0_real64
  real(real64), parameter :: max_aspect = 2000.0_real64
  real(real64), parameter :: min_area_ratio = 1.0e-6_real64

  real(real64), parameter :: te_wedge_rad = 0.3_real64   ! outer-ring half gap behind an airfoil trailing edge
  real(real64), parameter :: ilu_pivot_min = 1.0e-14_real64
  real(real64), parameter :: breakdown_min = 1.0e-30_real64
  real(real64), parameter :: kutta_den_min = 1.0e-12_real64
  real(real64), parameter :: cl_gap_flag = 0.05_real64
end module constants
