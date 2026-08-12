program test_force_and_ghosts
use CONFIG
use ALLDATA
use BOUNDARY
use FLUX
implicit none

type(cell_interface) :: face
real(8) :: fn,ft
integer :: i,j

phi=0.1d0
case_id=5
face%theta=0.37d0
face%cosa=1.d0
face%sina=0.d0
call face_acceleration(face,fn,ft)
call require(abs(fn-phi)<1.d-14 .and. abs(ft)<1.d-14, &
             'Poiseuille vertical face acceleration is not purely normal')
face%cosa=0.d0
face%sina=1.d0
call face_acceleration(face,fn,ft)
call require(abs(fn)<1.d-14 .and. abs(ft+phi)<1.d-14, &
             'Poiseuille horizontal face has a wall-normal acceleration')

phi=-1.5d0
case_id=6
face%theta=0.37d0
face%cosa=1.d0
face%sina=0.d0
call face_acceleration(face,fn,ft)
call require(abs(fn-phi*cos(face%theta))<1.d-14 .and. &
             abs(ft-phi*sin(face%theta))<1.d-14, &
             'RT radial acceleration projection is incorrect')

ixmin=1; ixmax=2; iymin=1; iymax=2
unum=2; vnum=2
allocate(ctr(0:3,0:3))
do j=0,3
  do i=0,3
    allocate(ctr(i,j)%h(unum,vnum),ctr(i,j)%b(unum,vnum))
    allocate(ctr(i,j)%sh(unum,vnum,2),ctr(i,j)%sb(unum,vnum,2))
    ctr(i,j)%w=0.d0
    ctr(i,j)%h=0.d0
    ctr(i,j)%b=0.d0
    ctr(i,j)%sh=0.d0
    ctr(i,j)%sb=0.d0
  end do
end do

do j=iymin,iymax
  do i=ixmin,ixmax
    ctr(i,j)%w=dble(10*i+j)
    ctr(i,j)%h=dble(10*i+j)
    ctr(i,j)%b=-dble(10*i+j)
  end do
end do

call bc_periodic_2d(Idirc)
call require(all(ctr(0,1)%w==ctr(2,1)%w), 'left periodic ghost is incorrect')
call require(all(ctr(3,2)%h==ctr(1,2)%h), 'right periodic ghost is incorrect')

write(*,*) 'PASS: case-aware force direction and periodic ghosts'

contains

subroutine require(condition,message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ',trim(message)
    error stop 1
  end if
end subroutine require

end program test_force_and_ghosts
