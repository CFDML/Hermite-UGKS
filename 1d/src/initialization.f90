module INITIALIZER
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
use VELOCITY
implicit none
contains
! ------------------------------------------------------------
	subroutine initialize()

	call allocate_alldata()
	call init_loopSetup()
	call init_geometry()
	!flow field initialization
	if(isRestart .eqv. .true.) then
		!call restart(fileInitId)
	else
		call init_flowField()
    end if

    call init_velocity()
	call init_distribution_function()

	end subroutine initialize
! ------------------------------------------------------------
	subroutine allocate_alldata()

	write(*,*) "allocating alldata"

	allocate(ctr(ixmin:ixmax)) !with ghost cell
	allocate(face(ixmax+1))

	!gravity field
	!重力定义在node point
	allocate(pot(1:ixmax+1))

	end subroutine allocate_alldata
! ------------------------------------------------------------
	subroutine init_loopSetup()

	! all run-time parameters (cfl, maxtime, knudsen, phi, TwL, TwR, case_id,
	! refpfile, ...) come from module CONFIG via input.namelist.
	interp_method=secondOrder
	simtime=0.d0
	iter=1

	! boundary type derived from case_id (not user-set):
	!   case 4 (Fourier) -> diffuse isothermal wall; cases 2,3 -> specular
	select case(case_id)
	case(2)
		bc_type=1
	case(3)
		bc_type=1
	case(4)
		bc_type=2
	case default
		write(*,*) 'Error: unknown case_id ',case_id
		stop
	end select

	end subroutine init_loopSetup
! ------------------------------------------------------------
	subroutine init_geometry
	integer:: i
	real(kind=8):: dx

	dx=xlength/(ixmax-ixmin+1)

	do i=ixmin,ixmax
		ctr(i)%x=(i-0.5d0)*dx
		ctr(i)%length=dx
	end do

	do i=1,ixmax+1
		face(i)%bcType=1
	end do
	write(*,*) "dx: ",dx

	!**********************************************
	!external force phi comes from CONFIG (input.namelist)
	do i=ixmin,ixmax+1
		pot(i)=-phi*dx*(i-1)
	end do

	end subroutine init_geometry
! ------------------------------------------------------------
	! ------------------------------------------------------------
	subroutine init_flowField()
	integer:: i,j

	write(*,*) "initializing flow field"
	gasR=R0/NAmass*1.d3
	gamma=get_gamma(innerK)

	rhoInf=1.d0
	velInf=0.d0
	! Tref comes from CONFIG (input.namelist); Tinf/lambdainf follow it
	Tinf=Tref
	lambdainf=1.d0/Tinf
	sosInf=1.d0

	momInf=rhoInf*velInf
	EInf=0.5d0*momInf*momInf/rhoInf+&       ! kinetic energy
		  0.5d0*rhoInf/lambdainf/(gamma-1.0)  ! internal energy

	omega=0.72d0
	alpha_ref=1.0
	omega_ref=0.5
	muref=get_mu(knudsen,alpha_ref,omega_ref) ! viscosity coefficient

	write(*,*) "gamma: ",gamma
	write(*,*) "refer viscosity: ",muref

	select case(case_id)
	!----------------------------------------------------------
	case(2)  ! 1D hydrostatic equilibrium with pressure perturbation
	!   rho0(x)=exp(ax), p0(x)=exp(ax), a=phi=-1, plus
	!   p = p0(x) + 0.01*exp(-100*(x-0.5)^2);  prim(3)=lambda=0.5*rho/p
	!----------------------------------------------------------
	do i=ixmin,ixmax
		ctr(i)%prim(1)=exp(-ctr(i)%x)
		ctr(i)%prim(2)=0.0d0
		ctr(i)%prim(3)=0.5d0*ctr(i)%prim(1)/&
			((exp(-ctr(i)%x))+0.01d0*exp(-100.d0*(ctr(i)%x-0.5d0)**2))
		ctr(i)%w=get_conserved(ctr(i)%prim)
		ctr(i)%sw=0.d0
	end do

	!----------------------------------------------------------
	case(3)  ! Sod shock tube: left (rho,U,p)=(1,0,1), right=(0.125,0,0.1)
	!   prim(3)=lambda=0.5*rho/p
	!----------------------------------------------------------
	do i=ixmin,ixmax/2
		ctr(i)%prim(1)=1.0d0
		ctr(i)%prim(2)=0.0d0
		ctr(i)%prim(3)=0.5d0          ! 0.5*1/1
		ctr(i)%w=get_conserved(ctr(i)%prim)
		ctr(i)%sw=0.d0
	end do
	do i=ixmax/2+1,ixmax
		ctr(i)%prim(1)=0.125d0
		ctr(i)%prim(2)=0.0d0
		ctr(i)%prim(3)=0.625d0        ! 0.5*0.125/0.1
		ctr(i)%w=get_conserved(ctr(i)%prim)
		ctr(i)%sw=0.d0
	end do

	!----------------------------------------------------------
	case(4)  ! Fourier flow: linear temperature ramp TwL->TwR, rho=1, U=0
	!   prim(3)=lambda=1/(2T)
	!----------------------------------------------------------
	do i=ixmin,ixmax
		ctr(i)%prim(1)=1.0d0
		ctr(i)%prim(2)=0.0d0
		ctr(i)%prim(3)=get_lambda_from_temperature(&
			TwL+(ctr(i)%x)*(TwR-TwL))
		ctr(i)%w=get_conserved(ctr(i)%prim)
		ctr(i)%sw=0.d0
	end do

	end select

	end subroutine init_flowField
! ------------------------------------------------------------
	subroutine init_distribution_function()
	integer:: i,j
	real(kind=8),allocatable,dimension(:):: H,B

	allocate(H(unum))
	allocate(B(unum))

	do i=ixmin,ixmax
		allocate(ctr(i)%h(unum))
		allocate(ctr(i)%b(unum))
		allocate(ctr(i)%sh(unum))
		allocate(ctr(i)%sb(unum))
	end do

	do i=1,ixmax+1
		allocate(face(i)%flux_h(unum))
		allocate(face(i)%flux_b(unum))
	end do

	do i=ixmin,ixmax
		call reduced_maxwell(H,B,ctr(i)%prim)
		ctr(i)%h=H
		ctr(i)%b=B
		ctr(i)%sh=0.d0
		ctr(i)%sb=0.d0
      !  print *, "H at index", i, "is:",H(1)
	end do

	end subroutine init_distribution_function
! ------------------------------------------------------------
end module INITIALIZER
