program test_csr_poisson
  use kinds
  use csr_solver
  use test_support
  implicit none

  type(csr_matrix) :: a, lu, empty, sing
  real(real64), allocatable :: xe(:), b(:), x(:), dense(:, :), yd(:), ym(:)
  real(real64) :: resid
  integer :: nx, n, it_pre, it_none, ierr, i, k
  integer, parameter :: nsmall = 4

  ! small case: matvec against dense multiply
  call poisson(nsmall, a)
  n = a%n
  allocate (dense(n, n), xe(n), yd(n), ym(n))
  dense = 0.0_real64
  do i = 1, n
    do k = a%rowptr(i), a%rowptr(i + 1) - 1
      dense(i, a%colidx(k)) = a%val(k)
    end do
    xe(i) = sin(0.7_real64 * i)
  end do
  yd = matmul(dense, xe)
  call csr_matvec(a, xe, ym)
  call check_close(maxval(abs(yd - ym)), 0.0_real64, 1.0e-14_real64, 'matvec vs dense')
  deallocate (dense, xe, yd, ym)

  ! solve with known solution
  nx = 40
  call poisson(nx, a)
  n = a%n
  allocate (xe(n), b(n), x(n), ym(n))
  do i = 1, n
    xe(i) = sin(0.01_real64 * i) + 0.3_real64
  end do
  call csr_matvec(a, xe, b)

  call ilu0_factor(a, lu, ierr)
  call check(ierr == ierr_ok, 'ilu0 factor ok')
  call bicgstab(a, lu, b, 1.0e-10_real64, 1000, .false., x, it_pre, resid, ierr)
  call check(ierr == ierr_ok, 'bicgstab ilu converges')
  call check(resid <= 1.0e-10_real64, 'true residual meets tol')
  call check_close(maxval(abs(x - xe)), 0.0_real64, 1.0e-7_real64, 'solution error')

  call bicgstab(a, empty, b, 1.0e-10_real64, 5000, .false., x, it_none, resid, ierr)
  call check(ierr == ierr_ok, 'bicgstab unpreconditioned converges')
  print '(a,i0,a,i0)', '  iterations: ilu0=', it_pre, ' none=', it_none
  call check(it_pre < it_none, 'ilu0 needs fewer iterations')

  ! stalled: iteration cap far too small
  call bicgstab(a, lu, b, 1.0e-14_real64, 2, .false., x, it_pre, resid, ierr)
  call check(ierr == ierr_noconv, 'iteration cap gives ierr_noconv')

  ! zero right-hand side
  b = 0.0_real64
  call bicgstab(a, lu, b, 1.0e-10_real64, 100, .false., x, it_pre, resid, ierr)
  call check(ierr == ierr_ok .and. all(x == 0.0_real64) .and. resid == 0.0_real64, 'zero rhs gives x=0')

  ! singular matrix [[1,1],[1,1]] and zero pivot [[0,1],[1,1]]
  call two_by_two(1.0_real64, sing)
  call ilu0_factor(sing, lu, ierr)
  call check(ierr == ierr_pivot, 'singular matrix gives ierr_pivot')
  call two_by_two(0.0_real64, sing)
  call ilu0_factor(sing, lu, ierr)
  call check(ierr == ierr_pivot, 'zero pivot gives ierr_pivot')

  ! breakdown: A = [[0,1],[-1,0]], b=(1,0): rhat.v = 0 on the first step
  call rotation(sing)
  b = 0.0_real64
  deallocate (x)
  allocate (x(2))
  call bicgstab(sing, empty, [1.0_real64, 0.0_real64], 1.0e-10_real64, 50, .false., x, it_pre, resid, ierr)
  call check(ierr == ierr_breakdown, 'rotation matrix gives ierr_breakdown')

  call finish()

contains

  subroutine poisson(m, mat)
    integer, intent(in) :: m
    type(csr_matrix), intent(out) :: mat
    integer :: ii, jj, row, kk
    mat%n = m * m
    allocate (mat%rowptr(mat%n + 1), mat%colidx(5 * mat%n), mat%val(5 * mat%n))
    kk = 0
    row = 0
    do jj = 1, m
      do ii = 1, m
        row = row + 1
        mat%rowptr(row) = kk + 1
        if (jj > 1) call put(mat, kk, row - m, -1.0_real64)
        if (ii > 1) call put(mat, kk, row - 1, -1.0_real64)
        call put(mat, kk, row, 4.0_real64)
        if (ii < m) call put(mat, kk, row + 1, -1.0_real64)
        if (jj < m) call put(mat, kk, row + m, -1.0_real64)
      end do
    end do
    mat%rowptr(mat%n + 1) = kk + 1
    mat%nnz = kk
  end subroutine poisson

  subroutine put(mat, kk, c, v)
    type(csr_matrix), intent(inout) :: mat
    integer, intent(inout) :: kk
    integer, intent(in) :: c
    real(real64), intent(in) :: v
    kk = kk + 1
    mat%colidx(kk) = c
    mat%val(kk) = v
  end subroutine put

  subroutine two_by_two(a11, mat)
    real(real64), intent(in) :: a11
    type(csr_matrix), intent(out) :: mat
    mat%n = 2
    mat%nnz = 4
    mat%rowptr = [1, 3, 5]
    mat%colidx = [1, 2, 1, 2]
    mat%val = [a11, 1.0_real64, 1.0_real64, 1.0_real64]
  end subroutine two_by_two

  subroutine rotation(mat)
    type(csr_matrix), intent(out) :: mat
    mat%n = 2
    mat%nnz = 2
    mat%rowptr = [1, 2, 3]
    mat%colidx = [2, 1]
    mat%val = [1.0_real64, -1.0_real64]
  end subroutine rotation

end program test_csr_poisson
