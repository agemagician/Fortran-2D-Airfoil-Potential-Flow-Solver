! Velocity, pressure coefficient and (slice 2) pressure forces.
module postproc
  use kinds, only: real64
  use constants, only: uinf
  use grid, only: grid_record, metrics_record
  implicit none
  private
  public :: velocity_field, surface_cp, cp_field, force_report, pressure_forces

  type :: force_report
    real(real64) :: cl_press = 0.0_real64, cd = 0.0_real64, cl_gap_rel = 0.0_real64
  end type force_report

contains

  ! u = psi_y, v = -psi_x from the chain rule; one-sided second-order eta differences at the rings.
  subroutine velocity_field(met, psi, u, v)
    type(metrics_record), intent(in) :: met
    real(real64), intent(in) :: psi(:, :)
    real(real64), intent(out) :: u(:, :), v(:, :)
    integer :: ni, nj, i, j, ip, im
    real(real64) :: p_xi, p_eta

    ni = size(psi, 1)
    nj = size(psi, 2)
    do j = 1, nj
      do i = 1, ni
        ip = mod(i, ni) + 1
        im = mod(i - 2 + ni, ni) + 1
        p_xi = 0.5_real64 * (psi(ip, j) - psi(im, j))
        if (j == 1) then
          p_eta = 0.5_real64 * (-3.0_real64 * psi(i, 1) + 4.0_real64 * psi(i, 2) - psi(i, 3))
        else if (j == nj) then
          p_eta = 0.5_real64 * (3.0_real64 * psi(i, nj) - 4.0_real64 * psi(i, nj - 1) + psi(i, nj - 2))
        else
          p_eta = 0.5_real64 * (psi(i, j + 1) - psi(i, j - 1))
        end if
        u(i, j) = (met%x_xi(i, j) * p_eta - met%x_eta(i, j) * p_xi) / met%jac(i, j)
        v(i, j) = (met%y_xi(i, j) * p_eta - met%y_eta(i, j) * p_xi) / met%jac(i, j)
      end do
    end do
  end subroutine velocity_field

  elemental function cp_field(u, v) result(cp)
    real(real64), intent(in) :: u, v
    real(real64) :: cp
    cp = 1.0_real64 - (u * u + v * v) / (uinf * uinf)
  end function cp_field

  ! Cp at the body-ring nodes (j = 1).
  subroutine surface_cp(u, v, cp)
    real(real64), intent(in) :: u(:, :), v(:, :)
    real(real64), intent(out) :: cp(:)
    cp = cp_field(u(:, 1), v(:, 1))
  end subroutine surface_cp

  ! Pressure force on the ring including the blunt face (segment ni -> 1). Cp_mid is the mean of the end nodes.
  ! cl_circ (CL from circulation) is only used for the relative gap.
  subroutine pressure_forces(g, cp_ring, alpha_rad, cl_circ, rep)
    type(grid_record), intent(in) :: g
    real(real64), intent(in) :: cp_ring(:), alpha_rad, cl_circ
    type(force_report), intent(out) :: rep
    real(real64) :: cx, cy, cpm
    integer :: i, ip

    cx = 0.0_real64
    cy = 0.0_real64
    do i = 1, g%ni
      ip = mod(i, g%ni) + 1
      cpm = 0.5_real64 * (cp_ring(i) + cp_ring(ip))
      cx = cx + cpm * (g%y(ip, 1) - g%y(i, 1))
      cy = cy - cpm * (g%x(ip, 1) - g%x(i, 1))
    end do
    rep%cl_press = cy * cos(alpha_rad) - cx * sin(alpha_rad)
    rep%cd = cx * cos(alpha_rad) + cy * sin(alpha_rad)
    rep%cl_gap_rel = abs(rep%cl_press - cl_circ) / max(abs(cl_circ), tiny(1.0_real64))
  end subroutine pressure_forces

end module postproc
