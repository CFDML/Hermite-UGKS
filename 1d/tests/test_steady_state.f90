program test_steady_state
use SOLVER, only: fourier_is_steady
implicit none

real(8) :: residual(3)

residual = [1.d-10, 2.d-10, 5.d-10]
call require(fourier_is_steady(residual, 2.d-12, 1.d-9, 1.d-10), &
             'converged state was rejected')

residual(2) = 2.d-8
call require(.not. fourier_is_steady(residual, 2.d-12, 1.d-9, 1.d-10), &
             'large macro residual was accepted')

residual = [1.d-10, 2.d-10, 5.d-10]
call require(.not. fourier_is_steady(residual, 2.d-7, 1.d-9, 1.d-10), &
             'large mass drift was accepted')

write(*,*) 'PASS: Fourier steady-state predicate'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_steady_state
