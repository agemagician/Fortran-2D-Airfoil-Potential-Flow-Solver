! Starting values, unmeasured until phase 7. Edit once from the measured study.
module tolerances
  use kinds, only: real64
  implicit none
  public
  real(real64), parameter :: tol_patch = 1.0e-12_real64
  real(real64), parameter :: tol_velocity_exact = 1.0e-10_real64
  real(real64), parameter :: tol_zero_cl = 1.0e-8_real64
  real(real64), parameter :: tol_gammaaux_rel = 1.0e-2_real64
  real(real64), parameter :: tol_cyl_cp_abs = 8.0e-2_real64
  real(real64), parameter :: tol_cyl_cl_rel = 1.0e-2_real64
  real(real64), parameter :: uniform_refine_ratio_min = 3.0_real64
  real(real64), parameter :: tol_gamma_change_rel = 0.02_real64
  real(real64), parameter :: joukowski_ratio_max = 0.95_real64   ! measured 0.89; initial estimate 0.8 assumed no Kutta bias
  real(real64), parameter :: tol_joukowski_extrap_rel = 0.15_real64   ! measured 0.119 (Kutta bias, see README); initial estimate 0.02
  real(real64), parameter :: tol_alpha_l0_deg = 0.5_real64
  real(real64), parameter :: tol_rotated_seam_cl = 1.0e-8_real64
  real(real64), parameter :: tol_naca_face_rel = 0.05_real64
end module tolerances
