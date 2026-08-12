module SOLVER
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
use FLUX
use GSOURCE


implicit none

real(kind=8) :: steady_res(DIM+2)=huge(1.d0)
    
contains
! ------------------------------------------------------------
	subroutine timestep()
    real(kind=8):: tmin !minimum dt allowed
    real(kind=8):: temp !max 1/dt allowed
    real(kind=8):: sos !speed of sound
    real(kind=8):: prim(DIM+2)
    integer:: i,j
    
    !set initial value
    tmin=0.0
    temp=0.0
    
    !$omp parallel
    !$omp do private(i,sos,prim) reduction(max:temp)
    do i=ixmin,ixmax
        prim=get_primary(ctr(i)%w)
        sos=get_sos(prim(3))

        !maximum velocity(plus speed of sound)
        prim(2)=max(umax,abs(prim(2)))+sos
        !prim(2)=abs(prim(2))+sos

        !maximum 1/dt allowed
        temp=max(temp,prim(2)/ctr(i)%length,abs(phi)/du)
    end do
    !$omp end do
    !$omp end parallel

    tmin=1.0/temp
    !time step
    dt=cfl*tmin/1.1d0
  !  print '(F10.5)', dt
	end subroutine timestep
! ------------------------------------------------------------
	subroutine interpolation()
    integer :: i

    call interp_boundary(ctr(ixmin),ctr(ixmin),ctr(ixmin+1))
    call interp_boundary(ctr(ixmax-1),ctr(ixmax),ctr(ixmax))
    !call interp_inner(ctr(ixmin),ctr(ixmin),ctr(ixmin+1))
    !call interp_inner(ctr(ixmax-1),ctr(ixmax),ctr(ixmax))

    !$omp parallel
    !$omp do
    do i=ixmin+1,ixmax-1
        call interp_inner(ctr(i-1),ctr(i),ctr(i+1))
    end do
    !omp end do nowait
    !$omp end parallel

    end subroutine interpolation
! ------------------------------------------------------------
	! one-sided interpolation of the boundary cell
	subroutine interp_boundary(cell_L,cell_N,cell_R)
    type(cell_center),intent(inout):: cell_N
    type(cell_center),intent(inout):: cell_L,cell_R

    cell_N%sh = (cell_R%h-cell_L%h)/(0.5*cell_R%length+0.5*cell_L%length)
    cell_N%sb = (cell_R%b-cell_L%b)/(0.5*cell_R%length+0.5*cell_L%length)

    !Conservative variable slope
    cell_N%sw=(cell_R%w-cell_L%w)/(0.5*cell_R%length+0.5*cell_L%length)
	end subroutine interp_boundary
! ------------------------------------------------------------
	! interpolation of the inner cell
	! using Van-Leer limiter
	subroutine interp_inner(cell_L,cell_N,cell_R)
    type(cell_center),intent(in):: cell_L,cell_R
    type(cell_center),intent(inout):: cell_N
    real(kind=8),allocatable,dimension(:):: sL,sR
    real(kind=8):: swL(DIM+2),swR(DIM+2)
    
    !allocate array
    allocate(sL(unum))
    allocate(sR(unum))

    sL=(cell_N%h-cell_L%h)/(0.5*cell_N%length+0.5*cell_L%length)
    sR=(cell_R%h-cell_N%h)/(0.5*cell_R%length+0.5*cell_N%length)
    cell_N%sh=(sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
    
    sL=(cell_N%b-cell_L%b)/(0.5*cell_N%length+0.5*cell_L%length)
    sR=(cell_R%b-cell_N%b)/(0.5*cell_R%length+0.5*cell_N%length)
    cell_N%sb=(sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
    
    !Conservative variable slope
    swL=(cell_N%w-cell_L%w)/(0.5*cell_N%length+0.5*cell_L%length)
    swR=(cell_R%w-cell_N%w)/(0.5*cell_R%length+0.5*cell_N%length)
    cell_N%sw=(sign(UP,swR)+sign(UP,swL))*abs(swR)*abs(swL)/(abs(swR)+abs(swL)+smv)

	end subroutine interp_inner
! ------------------------------------------------------------
	subroutine evolution()
	integer:: i
	! boundary condition at the two domain ends, selected by bc_type
	!   bc_type=1 : specular reflection (cases 2, 3)
	!   bc_type=2 : diffuse / isothermal Maxwellian wall (case 4, Fourier)
	! prim = (rho, u, lambda); lambda = 1/(2T). For a wall at temperature Tw,
	! rho is determined inside calc_flux_boundary by mass-flux balance (prim(1)=-R1/R2),
	! so the passed rho value is unused.
	real(kind=8) :: bcL(DIM+2), bcR(DIM+2)

	select case(bc_type)
	case(1)  ! specular reflection
		call calc_flux_specular(face(ixmin),ctr(ixmin),RN)
		call calc_flux_specular(face(ixmax+1),ctr(ixmax),RY)
	case(2)  ! diffuse isothermal wall
		bcL = (/0.d0, 0.d0, get_lambda_from_temperature(TwL)/)
		bcR = (/0.d0, 0.d0, get_lambda_from_temperature(TwR)/)
		call calc_flux_boundary(bcL,face(ixmin),ctr(ixmin),RN)
		call calc_flux_boundary(bcR,face(ixmax+1),ctr(ixmax),RY)
	case default
		write(*,*) "Error: unknown bc_type ",bc_type
		stop
	end select

    !$omp parallel
    !$omp do
    do i=ixmin+1,ixmax
        call calc_flux(ctr(i-1),face(i),ctr(i))
    end do

    !$omp end do nowait
    !$omp end parallel

	end subroutine evolution
! ------------------------------------------------------------
	subroutine update()
	!t=t^n
    real(kind=8),allocatable,dimension(:):: H_old,B_old !equilibrium distribution at t=t^n
    real(kind=8):: w_old(DIM+2)
    real(kind=8):: prim_old(DIM+2)
    real(kind=8):: tau_old
    
    !t=t^n+1
    real(kind=8),allocatable,dimension(:):: H,B !equilibrium distribution at t=t^{n+1}
    real(kind=8),allocatable,dimension(:):: H_plus,B_plus !Shakhov part
    
    real(kind=8):: prim(DIM+2) !primary variables at t^n and t^{n+1}
    real(kind=8):: qf
    real(kind=8):: tau !collision time and t^n and t^{n+1}
    real(kind=8):: sum_res(DIM+2),sum_avg(DIM+2)
    real(kind=8):: local_steady(DIM+2)

    
    integer:: i,j
    real(kind=8),allocatable,dimension(:):: hsource
    real(kind=8),allocatable,dimension(:):: bsource

    !allocate arrays
    allocate(H_old(unum))
    allocate(B_old(unum))
    allocate(H(unum))
    allocate(B(unum))
    allocate(H_plus(unum))
    allocate(B_plus(unum))
    allocate(hsource(unum))
    allocate(bsource(unum))

    !set initial value
    res=0.0
    steady_res=0.d0
    sum_res=0.0
    sum_avg=0.0

    !write(*,*) "update"
    do i=ixmin,ixmax
   ! do i=ixmin,ixmax
        !--------------------------------------------------
        ! store W^n and calculate H^n,B^n,\tau^n
        !--------------------------------------------------
        w_old=ctr(i)%w !store W^n
        
        prim_old=get_primary(w_old) !convert to primary variables
        call reduced_maxwell(H_old,B_old,prim_old) !calculate Maxwellian
        tau_old=get_tau(prim_old) !calculate collision time \tau^n

        !--------------------------------------------------
        ! update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
        !--------------------------------------------------
        ctr(i)%w=ctr(i)%w+(face(i)%flux-face(i+1)%flux)/ctr(i)%length

        prim=get_primary(ctr(i)%w)

        !projection
        prim(2)=prim(2)+phi*dt
        !ctr(i)%w(2)=ctr(i)%w(1)*prim(2)
        ctr(i)%w=get_conserved(prim)
        local_steady=primitive_change_residual(prim_old,prim,Tref)
        steady_res=max(steady_res,local_steady)

        call reduced_maxwell(H,B,prim)
        tau=get_tau(prim)

        !--------------------------------------------------
        ! record residual
        !--------------------------------------------------
        sum_res=sum_res+(w_old-ctr(i)%w)**2
        sum_avg=sum_avg+abs(ctr(i)%w)

        !--------------------------------------------------
        ! Shakhov part
        !--------------------------------------------------
        !heat flux at t=t^n
        qf=get_heat_flux(ctr(i)%h,ctr(i)%b,prim_old) 

        !h^+ = H+H^+ at t=t^n
        call shakhov_part(H_old,B_old,qf,prim_old,H_plus,B_plus) !H^+ and B^+
        H_old=H_old+H_plus !h^+
        B_old=B_old+B_plus !b^+

        !h^+ = H+H^+ at t=t^{n+1}
        call shakhov_part(H,B,qf,prim,H_plus,B_plus)
        H=H+H_plus
        B=B+B_plus

        !--------------------------------------------------
        ! body-force source term: finite-difference or Hermite spectral
        !--------------------------------------------------
        select case(trim(force_type))
        case('finite_difference')
            call cal_fv1(ctr(i))
        case('hermite')
            N_order1=Nm
            call cal_fv(ctr(i),prim_old)
        case default
            write(*,*) "unknown force_type '",trim(force_type),"'"
            stop
        end select

        hsource=dvh*dt*phi
        bsource=dvb*dt*phi

        !--------------------------------------------------
        ! update distribution function
        !--------------------------------------------------
      !  ctr(i)%h=(ctr(i)%h+(face(i)%flux_h-face(i+1)%flux_h)/ctr(i)%length+&
      !              0.5*dt*(H/tau+(H_old-ctr(i)%h)/tau_old)-hsource)/(1.0+0.5*dt/tau)
      !  ctr(i)%b=(ctr(i)%b+(face(i)%flux_b-face(i+1)%flux_b)/ctr(i)%length+&
      !              0.5*dt*(B/tau+(B_old-ctr(i)%b)/tau_old)-bsource)/(1.0+0.5*dt/tau)

        ctr(i)%h=(ctr(i)%h+(face(i)%flux_h-face(i+1)%flux_h)/ctr(i)%length+&
                    dt*(H/tau)-hsource)/(1.0+dt/tau)
        ctr(i)%b=(ctr(i)%b+(face(i)%flux_b-face(i+1)%flux_b)/ctr(i)%length+&
                    dt*(B/tau)-bsource)/(1.0+dt/tau)
    end do

    !final residual
    res=sqrt(ixmax*sum_res)/(sum_avg+SMV)

	end subroutine update
! ------------------------------------------------------------
    pure logical function fourier_is_steady(macro_residual, relative_mass_error, &
                                             residual_tolerance, mass_tolerance)
    real(kind=8),intent(in):: macro_residual(:)
    real(kind=8),intent(in):: relative_mass_error
    real(kind=8),intent(in):: residual_tolerance,mass_tolerance

    fourier_is_steady = maxval(abs(macro_residual)) <= residual_tolerance .and. &
                        abs(relative_mass_error) <= mass_tolerance
    end function fourier_is_steady
! ------------------------------------------------------------
    pure function primitive_change_residual(old_prim,new_prim,reference_temperature)
    real(kind=8),intent(in):: old_prim(DIM+2),new_prim(DIM+2)
    real(kind=8),intent(in):: reference_temperature
    real(kind=8):: primitive_change_residual(DIM+2)
    real(kind=8):: density_scale,velocity_scale,temperature_scale

    density_scale=max(1.d0,abs(old_prim(1)),abs(new_prim(1)))
    velocity_scale=dsqrt(reference_temperature)
    temperature_scale=reference_temperature

    primitive_change_residual(1)=abs(new_prim(1)-old_prim(1))/density_scale
    primitive_change_residual(2)=abs(new_prim(2)-old_prim(2))/velocity_scale
    primitive_change_residual(3)=abs(get_temperature_from_lambda(new_prim(3))-&
                                    get_temperature_from_lambda(old_prim(3)))/&
                                    temperature_scale
    end function primitive_change_residual
! ------------------------------------------------------------

end module SOLVER
