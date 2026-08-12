program test_gauss_hermite
use CONFIG
use ALLDATA
use VELOCITY
implicit none

integer, parameter :: orders(3) = [22, 34, 36]
integer :: k
real(8) :: symmetry_error, mass, second_moment

quad_type = 'gauss_hermite'
Tref = 1.d0

do k = 1, size(orders)
  nv = orders(k)
  if (allocated(uspace)) deallocate(uspace)
  if (allocated(weight)) deallocate(weight)

  call init_velocity()

  call require(all(uspace(2:) > uspace(:unum-1)), 'nodes are not strictly increasing')
  symmetry_error = maxval(abs(uspace + uspace(unum:1:-1)))
  call require(symmetry_error < 1.d-13, 'nodes are not symmetric')
  call require(all(weight > 0.d0), 'weights are not positive')

  mass = sum(weight*exp(-uspace**2/(2.d0*Tref))/sqrt(2.d0*PI*Tref))
  second_moment = sum(weight*uspace**2*exp(-uspace**2/(2.d0*Tref)) &
                      /sqrt(2.d0*PI*Tref))
  call require(abs(mass - 1.d0) < 5.d-10, 'Maxwell mass is inaccurate')
  call require(abs(second_moment - Tref) < 2.d-9, &
               'Maxwell second moment is inaccurate')
end do

write(*,*) 'PASS: Gauss-Hermite nodes, weights, and Maxwell moments'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message), ', Nv=', nv
    error stop 1
  end if
end subroutine require

end program test_gauss_hermite
