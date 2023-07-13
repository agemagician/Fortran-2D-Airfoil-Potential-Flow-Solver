module test_support
  use kinds, only: real64
  implicit none
  private
  public :: check, check_close, finish
  integer, save :: n_fail = 0, n_run = 0
contains
  subroutine check(ok, name)
    logical, intent(in) :: ok
    character(*), intent(in) :: name
    n_run = n_run + 1
    if (ok) then
      print '(a,a)', '  ok   ', name
    else
      n_fail = n_fail + 1
      print '(a,a)', '  FAIL ', name
    end if
  end subroutine check

  subroutine check_close(val, ref, tol, name)
    real(real64), intent(in) :: val, ref, tol
    character(*), intent(in) :: name
    real(real64) :: err
    err = abs(val - ref)
    n_run = n_run + 1
    if (err <= tol) then
      print '(a,a,a,es10.3,a,es9.2)', '  ok   ', name, '  err=', err, '  tol=', tol
    else
      n_fail = n_fail + 1
      print '(a,a,a,es10.3,a,es9.2)', '  FAIL ', name, '  err=', err, '  tol=', tol
    end if
  end subroutine check_close

  subroutine finish()
    print '(i0,a,i0,a)', n_run - n_fail, ' of ', n_run, ' checks passed'
    if (n_fail > 0) error stop 1
  end subroutine finish
end module test_support
