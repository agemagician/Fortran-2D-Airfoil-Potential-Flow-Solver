! Compressed sparse row matrix, ILU(0) and BiCGSTAB. No fluid-dynamics knowledge.
module csr_solver
  use kinds, only: real64, ierr_ok, ierr_pivot, ierr_breakdown, ierr_noconv
  use constants, only: ilu_pivot_min, breakdown_min
  implicit none
  private
  public :: csr_matrix, csr_matvec, ilu0_factor, bicgstab

  type :: csr_matrix
    integer :: n = 0, nnz = 0
    integer, allocatable :: rowptr(:), colidx(:), diag(:)  ! diag is set by ilu0_factor
    real(real64), allocatable :: val(:)
  end type csr_matrix

contains

  subroutine csr_matvec(a, x, y)
    type(csr_matrix), intent(in) :: a
    real(real64), intent(in) :: x(:)
    real(real64), intent(out) :: y(:)
    integer :: i, k
    do i = 1, a%n
      y(i) = 0.0_real64
      do k = a%rowptr(i), a%rowptr(i + 1) - 1
        y(i) = y(i) + a%val(k) * x(a%colidx(k))
      end do
    end do
  end subroutine csr_matvec

  ! Columns must be sorted ascending in every row.
  subroutine ilu0_factor(a, lu, ierr)
    type(csr_matrix), intent(in) :: a
    type(csr_matrix), intent(out) :: lu
    integer, intent(out) :: ierr
    integer, allocatable :: pos(:)
    integer :: i, k, j, kk, jj, c
    ierr = ierr_ok
    lu = a
    allocate (lu%diag(a%n), pos(a%n))
    lu%diag = 0
    pos = 0
    do i = 1, a%n
      do kk = lu%rowptr(i), lu%rowptr(i + 1) - 1
        pos(lu%colidx(kk)) = kk
        if (lu%colidx(kk) == i) lu%diag(i) = kk
      end do
      do kk = lu%rowptr(i), lu%rowptr(i + 1) - 1
        k = lu%colidx(kk)
        if (k >= i) exit
        lu%val(kk) = lu%val(kk) / lu%val(lu%diag(k))
        do jj = lu%diag(k) + 1, lu%rowptr(k + 1) - 1
          j = lu%colidx(jj)
          c = pos(j)
          if (c /= 0) lu%val(c) = lu%val(c) - lu%val(kk) * lu%val(jj)
        end do
      end do
      if (lu%diag(i) == 0) then
        ierr = ierr_pivot
        return
      end if
      if (abs(lu%val(lu%diag(i))) < ilu_pivot_min) then
        ierr = ierr_pivot
        return
      end if
      do kk = lu%rowptr(i), lu%rowptr(i + 1) - 1
        pos(lu%colidx(kk)) = 0
      end do
    end do
  end subroutine ilu0_factor

  ! z = (LU)^-1 r. An empty lu (n == 0) means no preconditioner.
  subroutine apply_precond(lu, r, z)
    type(csr_matrix), intent(in) :: lu
    real(real64), intent(in) :: r(:)
    real(real64), intent(out) :: z(:)
    integer :: i, k
    if (lu%n == 0) then
      z = r
      return
    end if
    do i = 1, lu%n
      z(i) = r(i)
      do k = lu%rowptr(i), lu%diag(i) - 1
        z(i) = z(i) - lu%val(k) * z(lu%colidx(k))
      end do
    end do
    do i = lu%n, 1, -1
      do k = lu%diag(i) + 1, lu%rowptr(i + 1) - 1
        z(i) = z(i) - lu%val(k) * z(lu%colidx(k))
      end do
      z(i) = z(i) / lu%val(lu%diag(i))
    end do
  end subroutine apply_precond

  ! Right-preconditioned BiCGSTAB, zero initial guess. resid is the true ||b-Ax||/||b||.
  subroutine bicgstab(a, lu, b, tol, maxiter, verbose, x, iters, resid, ierr)
    type(csr_matrix), intent(in) :: a, lu
    real(real64), intent(in) :: b(:), tol
    integer, intent(in) :: maxiter
    logical, intent(in) :: verbose
    real(real64), intent(out) :: x(:), resid
    integer, intent(out) :: iters, ierr
    real(real64), allocatable :: r(:), rhat(:), p(:), v(:), phat(:), s(:), shat(:), t(:)
    real(real64) :: bnorm, rho, rho_new, alpha, omega, beta, den, tt
    logical :: done

    ierr = ierr_ok
    iters = 0
    x = 0.0_real64
    resid = 0.0_real64
    bnorm = sqrt(dot_product(b, b))
    if (bnorm == 0.0_real64) return

    allocate (r(a%n), rhat(a%n), p(a%n), v(a%n), phat(a%n), s(a%n), shat(a%n), t(a%n))
    r = b
    rhat = r
    p = 0.0_real64
    v = 0.0_real64
    rho = 1.0_real64
    alpha = 1.0_real64
    omega = 1.0_real64
    done = .false.

    do while (iters < maxiter .and. .not. done)
      iters = iters + 1
      rho_new = dot_product(rhat, r)
      if (abs(rho_new) < breakdown_min) then
        ierr = ierr_breakdown
        exit
      end if
      beta = (rho_new / rho) * (alpha / omega)
      p = r + beta * (p - omega * v)
      call apply_precond(lu, p, phat)
      call csr_matvec(a, phat, v)
      den = dot_product(rhat, v)
      if (abs(den) < breakdown_min) then
        ierr = ierr_breakdown
        exit
      end if
      alpha = rho_new / den
      s = r - alpha * v
      if (sqrt(dot_product(s, s)) / bnorm < tol) then
        x = x + alpha * phat
        done = .true.
      else
        call apply_precond(lu, s, shat)
        call csr_matvec(a, shat, t)
        tt = dot_product(t, t)
        if (tt < breakdown_min) then
          ierr = ierr_breakdown
          exit
        end if
        omega = dot_product(t, s) / tt
        if (abs(omega) < breakdown_min) then
          ierr = ierr_breakdown
          exit
        end if
        x = x + alpha * phat + omega * shat
        r = s - omega * t
        rho = rho_new
        done = sqrt(dot_product(r, r)) / bnorm < tol
      end if
      if (verbose .and. mod(iters, 25) == 0) &
        print '(a,i6,a,es12.5)', '  bicgstab iter ', iters, '  resid ', sqrt(dot_product(r, r)) / bnorm
    end do

    call csr_matvec(a, x, t)
    resid = sqrt(dot_product(b - t, b - t)) / bnorm
    if (ierr == ierr_ok .and. (.not. done .or. resid > 10.0_real64 * tol)) ierr = ierr_noconv
  end subroutine bicgstab

end module csr_solver
