program test_force_moments
use CONFIG
use ALLDATA
use VELOCITY
use GASTHEORY
use GSOURCE
implicit none

type(cell_center) :: cell
real(8) :: prim(DIM+2)
real(8) :: mass_residual, momentum_residual, energy_residual

nv = 34
Nm = 6
Tref = 1.d0
quad_type = 'gauss_hermite'
gamma = get_gamma(innerK)
call init_velocity()

allocate(cell%h(unum), cell%b(unum))
prim = [1.d0, 0.3d0, 0.5d0]
call reduced_maxwell(cell%h, cell%b, prim)
N_order1 = Nm
call cal_fv(cell, prim)

mass_residual = sum(weight*dvh)
momentum_residual = sum(weight*uspace*dvh) + prim(1)
energy_residual = 0.5d0*sum(weight*(uspace**2*dvh + dvb)) &
                  + prim(1)*prim(2)

call require(abs(mass_residual) < 1.d-10, 'force mass moment is inconsistent')
call require(abs(momentum_residual) < 1.d-9, &
             'force momentum moment is inconsistent')
call require(abs(energy_residual) < 1.d-9, 'force energy moment is inconsistent')

write(*,*) 'PASS: Hermite force mass, momentum, and energy moments'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,'(A,3ES24.16)') 'force residuals: ', mass_residual, &
                            momentum_residual, energy_residual
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_force_moments
