! ------------------------------------------------------------
! Unified 1D UGKS + Hermite-spectral force-term solver
! Reference: Liu, Li, Zhang, Xiao, "An efficient solution algorithm for
!   force-driven continuum and rarefied flows", arXiv:2507.10021 (2025)
!
! This single source tree covers three 1D benchmarks by a case switch:
!   case_id = 2 : one-dimensional hydrostatic equilibrium (well-balanced)
!   case_id = 3 : Sod shock tube under external force
!   case_id = 4 : Fourier flow between two isothermal plates (diffuse wall)
! The algorithm (UGKS integral flux + Hermite spectral force term cal_fv)
! is identical across cases; only boundary and initial conditions differ.
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
integer :: i
integer :: steady_count
character(len=16) :: argstr
character(len=260) :: cfgfile
real(kind=8) :: initial_mass,current_mass,relative_mass_error
logical :: steady_reached

! ---- read all run parameters from a single namelist file ----
! Default config file is input.namelist in the run directory.
! Override with: ugks1d <configfile>
cfgfile = 'input.namelist'
argc = command_argument_count()
if (argc >= 1) then
    call get_command_argument(1, cfgfile)
end if
call read_config(cfgfile)

write(*,*) "UGKS + Hermite 1D solver, case_id =", case_id

! ------------------------------------------------------------
! Initialize parameter setup
! ------------------------------------------------------------
isRestart = .false.
iniType = 1 ! 1:generate 2:read
gridfile = '1Dperturb.neu'
inifile = ''

! ------------------------------------------------------------
! Initialize grid 
! ------------------------------------------------------------
 CALL SYSTEM_CLOCK(COUNT_RATE=rate)
   ! 获取程序开始时的时钟值
 CALL SYSTEM_CLOCK(COUNT=start)
call gridCheck()
call gridInit()

! ------------------------------------------------------------
! Initialize flow field
! ------------------------------------------------------------
call initialize()

initial_mass=0.d0
do i=ixmin,ixmax
    initial_mass=initial_mass+ctr(i)%w(1)*ctr(i)%length
end do
steady_reached=.false.
steady_count=0

! ------------------------------------------------------------
! Output initial field
! ------------------------------------------------------------
call write_macro(iter)

! ------------------------------------------------------------
! Main loop
! ------------------------------------------------------------
write(*,*) "begin main loop"
do while(.true.)
    !结束循环判据
    if(simtime >= maxtime) exit

    call timestep()
    call interpolation()
    call evolution()
    call update()

    !if(all(res<eps)) exit
    

    if(mod(iter,10)==0) then
        if(case_id==4) then
            write(*,"(A24,I8,6E15.7)") "iteration,time,dt,dprim:",iter,simtime,dt,steady_res
        else
    	    write(*,"(A24,I8,6E15.7)") "iteration,simtime,dt,res:",iter,simtime,dt,res
        end if
    end if

    if(mod(iter,1000)==0) then
        !call write_macro(iter)
    end if

    iter=iter+1
    simtime=simtime+dt

    if(case_id==4) then
        current_mass=0.d0
        do i=ixmin,ixmax
            current_mass=current_mass+ctr(i)%w(1)*ctr(i)%length
        end do
        relative_mass_error=abs(current_mass-initial_mass)/(abs(initial_mass)+SMV)
        if(fourier_is_steady(steady_res,relative_mass_error,steady_tol,mass_tol)) then
            steady_count=steady_count+1
        else
            steady_count=0
        end if
        if(steady_count>=steady_steps) then
            steady_reached=.true.
            write(*,'(A,I8,3ES15.7)') 'Fourier steady state reached: ',iter, &
                simtime,maxval(abs(steady_res)),relative_mass_error
            exit
        end if
    end if
end do

write(*,*) 'end of mainLoop'
write(*,*) 'iteration: ', iter
if(case_id==4 .and. .not.steady_reached) then
    write(*,*) 'WARNING: Fourier case reached maxtime before steady tolerance'
end if
call write_macro(iter)
CALL SYSTEM_CLOCK(COUNT=finish)
elapsed = REAL(finish - start) / REAL(rate)
PRINT *, "Elapsed time: ", elapsed, " seconds"
stop
end program main
