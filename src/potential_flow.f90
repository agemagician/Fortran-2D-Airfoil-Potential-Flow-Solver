! Stream-function Laplace solve on the O-grid, Kutta condition and circulation.
module potential_flow
  use kinds, only: real64, ierr_ok, ierr_kutta, ierr_bad_input
  use constants, only: uinf_c => uinf, kutta_den_min
  use grid, only: grid_record, metrics_record
  use csr_solver, only: csr_matrix, ilu0_factor, bicgstab
  implicit none
  private
  public :: solver_params, flow_result
  public :: stencil_weights, assemble, build_rhs, build_rhs_dirichlet, wall_q, body_flux, kutta, solve_flow

  type :: solver_params
    real(real64) :: alpha_rad = 0.0_real64
    real(real64) :: uinf = uinf_c
    real(real64) :: tol = 1.0e-10_real64
    integer :: maxiter = 2000
    logical :: use_kutta = .true.
    real(real64) :: c_prescribed = 0.0_real64
    integer :: kutta_i_lower = 2
    integer :: kutta_i_upper = 0          ! 0 means ni - 1
    logical :: verbose = .false.          ! passed to bicgstab (residual every 25 iterations)
  end type solver_params

  type :: flow_result
    real(real64), allocatable :: psi_stream(:, :), psi_aux(:, :), psi(:, :)
    real(real64) :: c = 0.0_real64, gamma_aux = 0.0_real64, gamma = 0.0_real64, gamma_stream_flux = 0.0_real64
    integer :: iters_stream = 0, iters_aux = 0
    real(real64) :: resid_stream = 0.0_real64, resid_aux = 0.0_real64
    integer :: ierr = ierr_ok
  end type flow_result

contains

  ! Nine-point weights w(di,dj) of -(d/dxi(A psi_xi - B psi_eta) + d/deta(C psi_eta - B psi_xi)) at node (i,j),
  ! with A = alpha_g/J, B = beta_g/J, C = gamma_g/J averaged onto faces. Used by both assembler and rhs.
  subroutine stencil_weights(met, i, j, w)
    type(metrics_record), intent(in) :: met
    integer, intent(in) :: i, j
    real(real64), intent(out) :: w(-1:1, -1:1)
    integer :: ni, ip, im
    real(real64) :: ae, aw, cn, cs, be, bw, bn, bs, q

    ni = size(met%jac, 1)
    ip = mod(i, ni) + 1
    im = mod(i - 2 + ni, ni) + 1
    ae = 0.5_real64 * (met%alpha_g(i, j) / met%jac(i, j) + met%alpha_g(ip, j) / met%jac(ip, j))
    aw = 0.5_real64 * (met%alpha_g(im, j) / met%jac(im, j) + met%alpha_g(i, j) / met%jac(i, j))
    cn = 0.5_real64 * (met%gamma_g(i, j) / met%jac(i, j) + met%gamma_g(i, j + 1) / met%jac(i, j + 1))
    cs = 0.5_real64 * (met%gamma_g(i, j - 1) / met%jac(i, j - 1) + met%gamma_g(i, j) / met%jac(i, j))
    be = 0.5_real64 * (met%beta_g(i, j) / met%jac(i, j) + met%beta_g(ip, j) / met%jac(ip, j))
    bw = 0.5_real64 * (met%beta_g(im, j) / met%jac(im, j) + met%beta_g(i, j) / met%jac(i, j))
    bn = 0.5_real64 * (met%beta_g(i, j) / met%jac(i, j) + met%beta_g(i, j + 1) / met%jac(i, j + 1))
    bs = 0.5_real64 * (met%beta_g(i, j - 1) / met%jac(i, j - 1) + met%beta_g(i, j) / met%jac(i, j))

    w = 0.0_real64
    ! east flux
    w(1, 0) = w(1, 0) + ae;  w(0, 0) = w(0, 0) - ae
    q = 0.25_real64 * be
    w(0, 1) = w(0, 1) - q;  w(0, -1) = w(0, -1) + q;  w(1, 1) = w(1, 1) - q;  w(1, -1) = w(1, -1) + q
    ! minus west flux
    w(0, 0) = w(0, 0) - aw;  w(-1, 0) = w(-1, 0) + aw
    q = 0.25_real64 * bw
    w(-1, 1) = w(-1, 1) + q;  w(-1, -1) = w(-1, -1) - q;  w(0, 1) = w(0, 1) + q;  w(0, -1) = w(0, -1) - q
    ! north flux
    w(0, 1) = w(0, 1) + cn;  w(0, 0) = w(0, 0) - cn
    q = 0.25_real64 * bn
    w(1, 0) = w(1, 0) - q;  w(-1, 0) = w(-1, 0) + q;  w(1, 1) = w(1, 1) - q;  w(-1, 1) = w(-1, 1) + q
    ! minus south flux
    w(0, 0) = w(0, 0) - cs;  w(0, -1) = w(0, -1) + cs
    q = 0.25_real64 * bs
    w(1, 0) = w(1, 0) + q;  w(-1, 0) = w(-1, 0) - q;  w(1, -1) = w(1, -1) + q;  w(-1, -1) = w(-1, -1) - q
    w = -w
  end subroutine stencil_weights

  ! Unknowns are nodes 2 <= j <= nj-1; row = (j-2)*ni + i. Columns sorted ascending in every row.
  subroutine assemble(g, met, a)
    type(grid_record), intent(in) :: g
    type(metrics_record), intent(in) :: met
    type(csr_matrix), intent(out) :: a
    real(real64) :: w(-1:1, -1:1)
    integer :: ni, nj, i, j, dj, m, di, iw, k, row, order(3)

    ni = g%ni
    nj = g%nj
    a%n = ni * (nj - 2)
    a%nnz = 9 * ni * (nj - 2) - 6 * ni
    allocate (a%rowptr(a%n + 1), a%colidx(a%nnz), a%val(a%nnz))
    k = 0
    row = 0
    do j = 2, nj - 1
      do i = 1, ni
        row = row + 1
        a%rowptr(row) = k + 1
        call stencil_weights(met, i, j, w)
        if (i == 1) then
          order = [0, 1, -1]
        else if (i == ni) then
          order = [1, -1, 0]
        else
          order = [-1, 0, 1]
        end if
        do dj = -1, 1
          if (j + dj < 2 .or. j + dj > nj - 1) cycle
          do m = 1, 3
            di = order(m)
            iw = mod(i + di - 1 + ni, ni) + 1
            k = k + 1
            a%colidx(k) = (j + dj - 2) * ni + iw
            a%val(k) = w(di, dj)
          end do
        end do
      end do
    end do
    a%rowptr(a%n + 1) = k + 1
  end subroutine assemble

  ! Right-hand side for Dirichlet values psi_in(1:ni) on j = 1 and psi_out(1:ni) on j = nj.
  subroutine build_rhs_dirichlet(g, met, psi_in, psi_out, b)
    type(grid_record), intent(in) :: g
    type(metrics_record), intent(in) :: met
    real(real64), intent(in) :: psi_in(:), psi_out(:)
    real(real64), intent(out) :: b(:)
    real(real64) :: w(-1:1, -1:1)
    integer :: ni, nj, i, j, di, iw, row

    ni = g%ni
    nj = g%nj
    b = 0.0_real64
    do j = 2, nj - 1
      if (j /= 2 .and. j /= nj - 1) cycle
      do i = 1, ni
        row = (j - 2) * ni + i
        call stencil_weights(met, i, j, w)
        do di = -1, 1
          iw = mod(i + di - 1 + ni, ni) + 1
          if (j == 2) b(row) = b(row) - w(di, -1) * psi_in(iw)
          if (j == nj - 1) b(row) = b(row) - w(di, 1) * psi_out(iw)
        end do
      end do
    end do
  end subroutine build_rhs_dirichlet

  subroutine build_rhs(g, met, par, b_stream, b_aux)
    type(grid_record), intent(in) :: g
    type(metrics_record), intent(in) :: met
    type(solver_params), intent(in) :: par
    real(real64), intent(out) :: b_stream(:), b_aux(:)
    real(real64), allocatable :: zero(:), outer(:)
    allocate (zero(g%ni), outer(g%ni))
    zero = 0.0_real64
    outer = par%uinf * (g%y(:, g%nj) * cos(par%alpha_rad) - g%x(:, g%nj) * sin(par%alpha_rad))
    call build_rhs_dirichlet(g, met, zero, outer, b_stream)
    outer = 1.0_real64
    call build_rhs_dirichlet(g, met, zero, outer, b_aux)
  end subroutine build_rhs

  ! Clockwise tangential speed on the body ring: q = psi_eta * l_xi / J, psi_eta one-sided second order.
  subroutine wall_q(psi, met, q)
    real(real64), intent(in) :: psi(:, :)
    type(metrics_record), intent(in) :: met
    real(real64), intent(out) :: q(:)
    integer :: i
    do i = 1, size(psi, 1)
      q(i) = 0.5_real64 * (-3.0_real64 * psi(i, 1) + 4.0_real64 * psi(i, 2) - psi(i, 3)) &
             * hypot(met%x_xi(i, 1), met%y_xi(i, 1)) / met%jac(i, 1)
    end do
  end subroutine wall_q

  ! Trapezoid line integral of q ds over all ni ring segments, blunt face (node ni to node 1) included.
  function body_flux(g, q) result(f)
    type(grid_record), intent(in) :: g
    real(real64), intent(in) :: q(:)
    real(real64) :: f
    integer :: i, ip
    f = 0.0_real64
    do i = 1, g%ni
      ip = mod(i, g%ni) + 1
      f = f + 0.5_real64 * (q(i) + q(ip)) * hypot(g%x(ip, 1) - g%x(i, 1), g%y(ip, 1) - g%y(i, 1))
    end do
  end function body_flux

  subroutine kutta(g, q_stream, q_aux, par, c, gamma_aux, gamma, ierr)
    type(grid_record), intent(in) :: g
    real(real64), intent(in) :: q_stream(:), q_aux(:)
    type(solver_params), intent(in) :: par
    real(real64), intent(out) :: c, gamma_aux, gamma
    integer, intent(out) :: ierr
    integer :: il, iu
    real(real64) :: den

    ierr = ierr_ok
    il = par%kutta_i_lower
    iu = par%kutta_i_upper
    if (iu <= 0) iu = g%ni - 1
    gamma_aux = body_flux(g, q_aux)
    if (par%use_kutta) then
      den = q_aux(il) + q_aux(iu)
      if (abs(den) < kutta_den_min) then
        ierr = ierr_kutta
        c = 0.0_real64
        gamma = 0.0_real64
        return
      end if
      c = -(q_stream(il) + q_stream(iu)) / den
    else
      c = par%c_prescribed
    end if
    ! total body flux; with a finite outer ring the stream part carries circulation too
    gamma = c * gamma_aux + body_flux(g, q_stream)
  end subroutine kutta

  subroutine solve_flow(g, met, par, res)
    type(grid_record), intent(in) :: g
    type(metrics_record), intent(in) :: met
    type(solver_params), intent(in) :: par
    type(flow_result), intent(out) :: res
    type(csr_matrix) :: a, lu
    real(real64), allocatable :: b_s(:), b_a(:), x(:), q_s(:), q_a(:)
    integer :: ni, nj, j

    ni = g%ni
    nj = g%nj
    call assemble(g, met, a)
    call ilu0_factor(a, lu, res%ierr)
    if (res%ierr /= ierr_ok) return
    allocate (b_s(a%n), b_a(a%n), x(a%n), q_s(ni), q_a(ni))
    call build_rhs(g, met, par, b_s, b_a)
    allocate (res%psi_stream(ni, nj), res%psi_aux(ni, nj), res%psi(ni, nj))

    call bicgstab(a, lu, b_s, par%tol, par%maxiter, par%verbose, x, res%iters_stream, res%resid_stream, res%ierr)
    if (res%ierr /= ierr_ok) return
    res%psi_stream(:, 1) = 0.0_real64
    res%psi_stream(:, nj) = par%uinf * (g%y(:, nj) * cos(par%alpha_rad) - g%x(:, nj) * sin(par%alpha_rad))
    do j = 2, nj - 1
      res%psi_stream(:, j) = x((j - 2) * ni + 1:(j - 1) * ni)
    end do

    call bicgstab(a, lu, b_a, par%tol, par%maxiter, par%verbose, x, res%iters_aux, res%resid_aux, res%ierr)
    if (res%ierr /= ierr_ok) return
    res%psi_aux(:, 1) = 0.0_real64
    res%psi_aux(:, nj) = 1.0_real64
    do j = 2, nj - 1
      res%psi_aux(:, j) = x((j - 2) * ni + 1:(j - 1) * ni)
    end do

    call wall_q(res%psi_stream, met, q_s)
    call wall_q(res%psi_aux, met, q_a)
    res%gamma_stream_flux = body_flux(g, q_s)
    call kutta(g, q_s, q_a, par, res%c, res%gamma_aux, res%gamma, res%ierr)
    if (res%ierr /= ierr_ok) return
    res%psi = res%psi_stream + res%c * res%psi_aux
  end subroutine solve_flow

end module potential_flow
