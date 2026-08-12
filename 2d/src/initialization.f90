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
	integer:: i

	write(*,*) "allocating alldata"
	if(grid%numElem<=0) then
		write(*,*) "minus grid element number: ",grid%numElem
		stop
	end if

	! allocate with one ghost-cell layer on each side so periodic and
	! diffuse boundary handlers can fill ctr(ixmin-1), ctr(ixmax+1), etc.
	allocate(ctr(ixmin-1:ixmax+1,iymin-1:iymax+1))
	allocate(vface(ixmin:ixmax+1,iymin:iymax))
	allocate(hface(ixmin:ixmax,iymin:iymax+1))
	
	end subroutine allocate_alldata
! ------------------------------------------------------------
	subroutine init_loopSetup()

	! all run-time parameters (cfl, maxtime, knudsen, phi, TwL, TwR, case_id)
	! come from module CONFIG via input.namelist.
	interp_method=secondOrder
	simtime=0.d0
	time1=0.d0
	if (case_id==6 .and. maxtime<=0.2d0) then
		time2=0.04d0
		time3=0.08d0
		time4=0.14d0
	else
		time2=0.8d0
		time3=1.4d0
		time4=2.d0
	end if
	iter=1

	end subroutine init_loopSetup
! ------------------------------------------------------------
	subroutine init_geometry
	integer:: i,j
	real(kind=8):: dx,dy

	dx=xlength/(grid%xindex-1)
	dy=ylength/(grid%yindex-1)

	! set geometry for physical cells and ghost cells (ghost cells get
	! the same dx/dy/area; their x/y are extrapolated, used only for flux)
	do i=ixmin-1,ixmax+1
	do j=iymin-1,iymax+1
		ctr(i,j)%x=(i-0.5d0)*dx
		ctr(i,j)%y=(j-0.5d0)*dy
		ctr(i,j)%r=dsqrt(ctr(i,j)%x**2+ctr(i,j)%y**2)
		ctr(i,j)%theta=datan2(ctr(i,j)%y,ctr(i,j)%x)
		ctr(i,j)%area=dx*dy
		ctr(i,j)%length(1)=dx
		ctr(i,j)%length(2)=dy
	end do
	end do

	do i=ixmin,ixmax
	do j=iymin,iymax
		vface(i,j)%length=dy
		vface(i,j)%cosa=1.d0
		vface(i,j)%sina=0.d0
		vface(i,j)%theta=datan2(ctr(i,j)%y,ctr(i,j)%x-0.5*ctr(i,j)%length(1))
	end do
	end do
	do j=iymin,iymax
		vface(grid%xindex,j)%length=dy
		vface(grid%xindex,j)%cosa=1.d0
		vface(grid%xindex,j)%sina=0.d0
		vface(grid%xindex,j)%theta=datan2(ctr(ixmax,j)%y,ctr(ixmax,j)%x+0.5*ctr(ixmax,j)%length(1))
	end do

	do i=ixmin,ixmax
	do j=iymin,iymax
		hface(i,j)%length=dx
		hface(i,j)%cosa=0.d0
		hface(i,j)%sina=1.d0
		hface(i,j)%theta=datan2(ctr(i,j)%y-0.5*ctr(i,j)%length(2),ctr(i,j)%x)
	end do
	end do
	do i=ixmin,ixmax
		hface(i,grid%yindex)%length=dx
		hface(i,grid%yindex)%cosa=0.d0
		hface(i,grid%yindex)%sina=1.d0
		hface(i,grid%yindex)%theta=datan2(ctr(i,iymax)%y+0.5*ctr(i,iymax)%length(2),ctr(i,iymax)%x)
	end do

	end subroutine init_geometry
! ------------------------------------------------------------
	subroutine init_flowField()
	integer:: i,j,k

	write(*,*) "initializing flow field"
	gasR=R0/NAmass*1.d3
	gamma=get_gamma(innerK)

	! Tref comes from CONFIG (input.namelist); Tinf/lambdainf follow it
	Tinf=Tref
	lambdainf=1.d0/Tinf
	sosInf=1.d0

	omega=0.72d0
	alpha_ref=1.0
	omega_ref=0.5
	if (case_id==5) then
		! Poiseuille reference used for Fig. 8 in the supplied source.
		muref=0.5d0*sqrt(PI)*knudsen
	else
		! Rayleigh--Taylor hard-sphere/VHS normalization.
		muref=get_mu(knudsen,alpha_ref,omega_ref)
	end if

	write(*,*) "gas gamma: ",gamma
	write(*,*) "refer viscosity: ",muref

	select case(case_id)
	!----------------------------------------------------------
	case(6)  ! Rayleigh-Taylor: polar initial condition on Cartesian grid
	prandtl=2.d0/3.d0
	!   rho0(r)=exp(-alpha*(r+r0)), alpha=2.68 (r<=r1) / 5.53 (r>r1)
	!   density interface r1=0.6*(1+0.02*cos(20*theta))
	!   pressure interface r1=0.62324965; prim(4)=lambda=0.5*rho/p
	!----------------------------------------------------------
	do i=ixmin,ixmax
	do j=iymin,iymax
		radius=0.6d0*(1.d0+0.02d0*cos(20*ctr(i,j)%theta))
		if(ctr(i,j)%r<=radius) then
			ctr(i,j)%prim(1)=exp(-2.68d0*(ctr(i,j)%r+0.258d0))
		else
			ctr(i,j)%prim(1)=exp(-5.53d0*(ctr(i,j)%r-0.308d0))
		end if
		radius=0.62324965d0
		if(ctr(i,j)%r<=radius) then
			ctr(i,j)%prim(4)=0.5d0*ctr(i,j)%prim(1)/1.5d0/exp(-2.68d0*(ctr(i,j)%r+0.258d0))*2.68d0
		else
			ctr(i,j)%prim(4)=0.5d0*ctr(i,j)%prim(1)/1.5d0/exp(-5.53d0*(ctr(i,j)%r-0.308d0))*5.53d0
		end if
		ctr(i,j)%prim(2:3)=0.d0
		ctr(i,j)%w=get_conserved(ctr(i,j)%prim)
		ctr(i,j)%sw=0.d0
	end do
	end do

	!----------------------------------------------------------
	case(5)  ! Poiseuille: uniform stationary gas, force phi in +x
	prandtl=1.d0
	!   rho=1, U=V=0, T=Tref; Rc=0.5 and prim(4)=lambda=1/T
	!----------------------------------------------------------
	do i=ixmin,ixmax
	do j=iymin,iymax
		ctr(i,j)%prim(1)=1.d0
		ctr(i,j)%prim(2)=0.d0
		ctr(i,j)%prim(3)=0.d0
		ctr(i,j)%prim(4)=1.d0/Tref
		ctr(i,j)%w=get_conserved(ctr(i,j)%prim)
		ctr(i,j)%sw=0.d0
	end do
	end do

	case default
		write(*,*) "Error: unknown case_id ",case_id
		stop
	end select

	end subroutine init_flowField
! ------------------------------------------------------------
	subroutine init_distribution_function()
	integer:: i,j,k
	real(kind=8),allocatable,dimension(:,:):: H,B

	allocate(H(unum,vnum))
	allocate(B(unum,vnum))

	! allocate distribution arrays for physical cells AND ghost cells
	! (ghost cells needed by periodic / diffuse boundary handlers)
	do i=ixmin-1,ixmax+1
	do j=iymin-1,iymax+1
		allocate(ctr(i,j)%h(unum,vnum))
		allocate(ctr(i,j)%b(unum,vnum))
		allocate(ctr(i,j)%sh(unum,vnum,2))
		allocate(ctr(i,j)%sb(unum,vnum,2))
	end do
	end do

	do i=ixmin,ixmax+1
	do j=iymin,iymax
		allocate(vface(i,j)%flux_h(unum,vnum))
		allocate(vface(i,j)%flux_b(unum,vnum))
	end do
	end do

	do i=ixmin,ixmax
	do j=iymin,iymax+1
		allocate(hface(i,j)%flux_h(unum,vnum))
		allocate(hface(i,j)%flux_b(unum,vnum))
	end do
	end do

	do i=ixmin,ixmax
	do j=iymin,iymax
		call reduced_maxwell(H,B,uspace,vspace,ctr(i,j)%prim)
		ctr(i,j)%h=H
		ctr(i,j)%b=B
		ctr(i,j)%sh=0.d0
		ctr(i,j)%sb=0.d0
	end do
    end do
    
	end subroutine init_distribution_function
! ------------------------------------------------------------
end module INITIALIZER
