program test_diffuse_wall
use CONFIG
use ALLDATA
use VELOCITY
use GASTHEORY
use FLUX
implicit none

type(cell_center) :: cell
type(cell_interface) :: wall
real(8) :: gas_prim(DIM+2), wall_prim(DIM+2)
real(8) :: left_mass_flux, right_mass_flux

nv = 22
Tref = 1.d0
quad_type = 'gauss_hermite'
call init_velocity()
dt = 1.d0

allocate(cell%h(unum), cell%b(unum), cell%sh(unum), cell%sb(unum))
allocate(wall%flux_h(unum), wall%flux_b(unum))
cell%length = 1.d0
cell%sh = 0.d0
cell%sb = 0.d0
gas_prim = [1.d0, 0.d0, 0.5d0]
wall_prim = [0.d0, 0.d0, 0.5d0]
call reduced_maxwell(cell%h, cell%b, gas_prim)

call calc_flux_boundary(wall_prim, wall, cell, RN)
left_mass_flux = wall%flux(1)
call calc_flux_boundary(wall_prim, wall, cell, RY)
right_mass_flux = wall%flux(1)

call require(abs(left_mass_flux) < 1.d-12, 'left diffuse wall leaks mass')
call require(abs(right_mass_flux) < 1.d-12, 'right diffuse wall leaks mass')
call require(abs(left_mass_flux + right_mass_flux) < 1.d-12, &
             'paired diffuse walls do not balance mass')

write(*,*) 'PASS: left and right diffuse walls have zero discrete mass flux'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,'(A,2ES24.16)') 'wall mass fluxes: ', left_mass_flux, right_mass_flux
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_diffuse_wall
