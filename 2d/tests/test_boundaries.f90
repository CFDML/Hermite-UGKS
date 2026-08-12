program test_boundaries
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
use VELOCITY
use BOUNDARY
use FLUX
implicit none

type(cell_center) :: cell
type(cell_interface) :: face
real(8) :: prim(4), recovered
integer :: i,j

nv = 16
Tref = 1.d0
quad_type = 'gauss_hermite'
call init_velocity()

allocate(cell%h(unum,vnum),cell%b(unum,vnum))
allocate(cell%sh(unum,vnum,2),cell%sb(unum,vnum,2))
allocate(face%flux_h(unum,vnum),face%flux_b(unum,vnum))
cell%length = 1.d0
cell%sh = 0.d0
cell%sb = 0.d0
prim = [1.d0,0.17d0,-0.11d0,1.d0]
call reduced_maxwell(cell%h,cell%b,uspace,vspace,prim)
cell%h = cell%h*(1.d0+0.04d0*uspace+0.03d0*vspace)
cell%b = cell%b*(1.d0-0.02d0*uspace+0.01d0*vspace)

dt = 1.d0
face%length = 1.d0
face%cosa = 1.d0
face%sina = 0.d0
call calc_flux_specular_vertical(face,cell,Idirc,RN)
call require(abs(face%flux(1)) < 2.d-11, 'left symmetry has nonzero mass flux')
do j=1,vnum
  do i=1,unum
    recovered = face%flux_h(i,j)/uspace(i,j)
    if (uspace(i,j) < 0.d0) then
      call require(abs(recovered-cell%h(i,j)) < 2.d-12, &
                   'left symmetry changed outgoing half range')
    else
      call require(abs(recovered-cell%h(unum-i+1,j)) < 2.d-12, &
                   'left symmetry reflected incoming half range incorrectly')
    end if
  end do
end do
call calc_flux_specular_vertical(face,cell,Idirc,RY)
call require(abs(face%flux(1)) < 2.d-11, 'right specular wall has nonzero mass flux')

face%cosa = 0.d0
face%sina = 1.d0
call calc_flux_specular_horizon(face,cell,Jdirc,RN)
call require(abs(face%flux(1)) < 2.d-11, 'bottom symmetry has nonzero mass flux')
do j=1,vnum
  do i=1,unum
    recovered = face%flux_h(i,j)/vspace(i,j)
    if (vspace(i,j) < 0.d0) then
      call require(abs(recovered-cell%h(i,j)) < 2.d-12, &
                   'bottom symmetry changed outgoing half range')
    else
      call require(abs(recovered-cell%h(i,vnum-j+1)) < 2.d-12, &
                   'bottom symmetry reflected incoming half range incorrectly')
    end if
  end do
end do
call calc_flux_specular_horizon(face,cell,Jdirc,RY)
call require(abs(face%flux(1)) < 2.d-11, 'top specular wall has nonzero mass flux')

face%cosa = 0.d0
face%sina = 1.d0
call calc_flux_boundary([0.d0,0.d0,0.d0,1.d0],face,cell,Jdirc,RN)
call require(abs(face%flux(1)) < 2.d-11, 'diffuse wall has nonzero mass flux')
call calc_flux_boundary([0.d0,0.d0,0.d0,1.d0],face,cell,Jdirc,RY)
call require(abs(face%flux(1)) < 2.d-11, 'upper diffuse wall has nonzero mass flux')

write(*,*) 'PASS: diffuse and symmetry boundary half-range balances'

contains

subroutine require(condition,message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_boundaries
