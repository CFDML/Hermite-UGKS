program test_primitive_change
use SOLVER, only: primitive_change_residual
implicit none

real(8) :: old_prim(3), new_prim(3), change(3)

old_prim = [1.d0, 0.2d0, 0.5d0]
new_prim = [1.01d0, 0.21d0, 0.4d0]
change = primitive_change_residual(old_prim, new_prim, 2.d0)

call require(abs(change(1) - 0.01d0/1.01d0) < 1.d-14, &
             'density change has the wrong scale')
call require(abs(change(2) - 0.01d0/sqrt(2.d0)) < 1.d-14, &
             'velocity change has the wrong thermal scale')
call require(abs(change(3) - 0.25d0) < 1.d-14, &
             'temperature change has the wrong reference scale')

write(*,*) 'PASS: primitive steady changes use fixed physical scales'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_primitive_change
