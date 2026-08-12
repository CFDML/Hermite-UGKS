! ------------------------------------------------------------
! Unified 2D UGKS + Hermite-spectral force-term solver
! Reference: Liu, Li, Zhang, Xiao, arXiv:2507.10021 (2025)
!
! This single source tree covers two 2D benchmarks by a case switch:
!   case_id = 6 : Rayleigh-Taylor instability (left/bottom symmetry,
!                 right/top closed specular faces)
!   case_id = 5 : Poiseuille flow (x-periodic, y-diffuse walls)
! The algorithm (UGKS integral flux + Hermite spectral force term) is
! identical across cases; only boundary and initial conditions differ.
! Based on the UGKS-Gravity code of Tianbai Xiao (2016).
! ------------------------------------------------------------
program main
use CONFIG
use ALLDATA
use GRIDMAKER
use INITIALIZER
use SOLVER
use OUTPUT
implicit none
INTEGER :: start, finish, rate
REAL :: elapsed
integer :: argc
integer :: i,j
character(len=64) :: cfgfile
real(8) :: mass_initial,mass_final,max_u_axis,max_v_axis
real(8) :: prim_diag(DIM+2)
logical :: wrote_time2,wrote_time3,wrote_time4

! ---- read all run parameters from a single namelist file ----
cfgfile = 'input.namelist'
argc = command_argument_count()
if (argc >= 1) call get_command_argument(1, cfgfile)
call read_config(cfgfile)
write(*,*) "UGKS + Hermite 2D solver, case_id =", case_id

! ------------------------------------------------------------
! Initialize parameter setup
! ------------------------------------------------------------
isRestart = .false.
iniType = 1 ! 1:generate 2:read
gridfile = 'rtmesh.neu'
inifile = ''

! ------------------------------------------------------------
! Initialize grid 
! ------------------------------------------------------------
CALL SYSTEM_CLOCK(COUNT_RATE=rate)
   ! ��ȡ����ʼʱ��ʱ��ֵ
CALL SYSTEM_CLOCK(COUNT=start)
call gridCheck()
call gridInit()

! ------------------------------------------------------------
! Initialize flow field
! ------------------------------------------------------------
call initialize()

! ------------------------------------------------------------
! Output initial field
! ------------------------------------------------------------
call write_macro(iter)
mass_initial=0.d0
do j=iymin,iymax
  do i=ixmin,ixmax
    mass_initial=mass_initial+ctr(i,j)%w(1)*ctr(i,j)%area
  end do
end do
wrote_time2=.false.
wrote_time3=.false.
wrote_time4=.false.

! ------------------------------------------------------------
! Main loop
! ------------------------------------------------------------
write(*,*) "begin main loop"
do while(.true.)
    call timestep()
    call interpolation()
    call evolution()
    call update()
    simtime=simtime+dt

    if(mod(iter,10)==0) then
    	write(*,"(A24,I8,6E15.7)") "iteration,simtime,dt,res:",iter,simtime,dt,res
    end if

    if (case_id==6 .and. .not.wrote_time2 .and. simtime>=time2) then
      call write_macro(iter)
      wrote_time2=.true.
    end if
    if (case_id==6 .and. .not.wrote_time3 .and. simtime>=time3) then
      call write_macro(iter)
      wrote_time3=.true.
    end if
    if (case_id==6 .and. .not.wrote_time4 .and. simtime>=time4) then
      call write_macro(iter)
      wrote_time4=.true.
    end if
    if (case_id==5 .and. mod(iter,5000)==0) call write_macro(iter)

    iter=iter+1
    if(simtime >= maxtime) exit
    if(case_id==5 .and. all(res<eps)) exit
    !if(ctr(ixmax*55/100,iymin)%w(1)>=ctr(ixmax*65/100,iymin)%w(1)) then
    !stop
    !end if
end do

write(*,*) 'end of mainLoop'
write(*,*) 'iteration: ', iter
write(*,*) "equi time vs Knudsen: ",simtime,Knudsen
call write_macro(iter)
mass_final=0.d0
max_u_axis=0.d0
max_v_axis=0.d0
do j=iymin,iymax
  do i=ixmin,ixmax
    mass_final=mass_final+ctr(i,j)%w(1)*ctr(i,j)%area
  end do
  prim_diag=get_primary(ctr(ixmin,j)%w)
  max_u_axis=max(max_u_axis,abs(prim_diag(2)))
end do
do i=ixmin,ixmax
  prim_diag=get_primary(ctr(i,iymin)%w)
  max_v_axis=max(max_v_axis,abs(prim_diag(3)))
end do
write(*,'(A,ES24.16)') 'initial_mass = ',mass_initial
write(*,'(A,ES24.16)') 'final_mass   = ',mass_final
write(*,'(A,ES24.16)') 'relative_mass_change = ',(mass_final-mass_initial)/mass_initial
write(*,'(A,ES24.16)') 'max_abs_u_left_axis  = ',max_u_axis
write(*,'(A,ES24.16)') 'max_abs_v_bottom_axis= ',max_v_axis
CALL SYSTEM_CLOCK(COUNT=finish)
elapsed = REAL(finish - start) / REAL(rate)
PRINT *, "Elapsed time: ", elapsed, " seconds"
stop
end program main
