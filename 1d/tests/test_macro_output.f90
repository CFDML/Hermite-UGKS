program test_macro_output
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
use GSOURCE
use VELOCITY
use OUTPUT
implicit none

integer, parameter :: output_index = 900001
integer :: unit_id, i
character(len=256) :: header, zone
character(len=40) :: output_name
real(8) :: values(10)
real(8) :: prim(DIM+2)

nx = 1
ixmin = 1
ixmax = 1
case_id = 3
nv = 22
Tref = 1.d0
quad_type = 'gauss_hermite'
gamma = get_gamma(innerK)
error1 = 0.d0
error2 = 0.d0
N_order1 = 0

call init_velocity()
allocate(ctr(1))
allocate(ctr(1)%h(unum), ctr(1)%b(unum))
ctr(1)%x = 0.5d0
ctr(1)%length = 1.d0
prim = [1.d0, 0.3d0, 0.5d0]
ctr(1)%w = get_conserved(prim)
call reduced_maxwell(ctr(1)%h, ctr(1)%b, prim)

call write_macro(output_index)
write(output_name,'(A,I0,A)') 'macro', output_index, '.dat'
open(newunit=unit_id, file=trim(output_name), status='old', action='read')
read(unit_id,'(A)') header
call require(index(header, 'Q') > 0, 'macro header has no heat-flux column')
read(unit_id,'(A)') zone
do i = 1, size(values)
  read(unit_id,*) values(i)
end do
close(unit_id, status='delete')

call require(abs(values(4) - 2.d0) < 1.d-12, &
             'temperature is not T=1/lambda')
call require(abs(values(5) - 1.d0) < 1.d-12, &
             'pressure is not rho*T/2')
call require(abs(values(7)) < 3.d-8, &
             'equilibrium Maxwellian has nonzero heat flux')

write(*,*) 'PASS: macro output contains physical T, p, and heat flux'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,'(A,10ES14.5)') 'output values: ', values
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_macro_output
