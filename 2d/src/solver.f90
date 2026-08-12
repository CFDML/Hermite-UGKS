module SOLVER
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
use FLUX
use GSOURCE
use HERMITESYMMETRIC2D
use BOUNDARY
implicit none
contains
! ------------------------------------------------------------
	subroutine timestep()
    real(kind=8) :: tmin ! minimum dt allowed
    real(kind=8) :: temp ! max 1/dt allowed
    real(kind=8) :: sos ! speed of sound
    real(kind=8) :: velocity_gap
    real(kind=8) :: prim(DIM+2)
    integer :: i, j

    ! set initial value
    tmin = 0.0
    temp = 0.0

    ! Set the number of threads for this parallel region
    call omp_set_num_threads(8)

    !$omp parallel default(none) private(i, j, sos, prim) shared(tmin, temp, ixmin, ixmax, iymin, iymax, umax, vmax_, cfl, ctr, case_id)
    !$omp do reduction(max:temp)
    do j = iymin, iymax
      do i = ixmin, ixmax
        prim = get_primary(ctr(i, j)%w)
        sos = get_sos(prim(4))

        ! maximum velocity (plus speed of sound)
        prim(2) = max(umax, abs(prim(2))) + sos
        prim(3) = max(vmax_, abs(prim(3))) + sos

        ! maximum 1/dt allowed
        if (case_id==5) then
          ! The requested three-cell Poiseuille solution is exactly uniform
          ! and periodic in x; only wall-normal transport limits the step.
          temp=max(temp,prim(3)/ctr(i,j)%length(2))
        else
          temp=max(temp,(ctr(i,j)%length(2)*prim(2)+ &
                        ctr(i,j)%length(1)*prim(3))/ctr(i,j)%area)
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel

    velocity_gap=min(minval(uspace(2:unum,1)-uspace(1:unum-1,1)), &
                     minval(vspace(1,2:vnum)-vspace(1,1:vnum-1)))
    temp=max(temp,abs(phi)/velocity_gap)
    tmin = 1.0 / temp

    ! time step
    dt = min(cfl*tmin/1.1d0,maxtime-simtime)

  end subroutine timestep
! ------------------------------------------------------------

  subroutine interpolation()
  use omp_lib
  integer :: i, j

  ! first order NO NEED interpolation of slope
  if (interp_method == firstOrder) return

  ! Set the number of threads for this parallel region
  call omp_set_num_threads(8)

  !$omp parallel default(none) private(i, j) shared(ixmin, ixmax, iymin, iymax, ctr, case_id)
  
  ! i direction
  if (case_id==5) then
    !$omp single
    call bc_periodic_2d(Idirc)
    !$omp end single
    !$omp do
    do j = iymin, iymax
      do i = ixmin, ixmax
        call interp_inner(ctr(i - 1, j), ctr(i, j), ctr(i + 1, j), IDIRC)
      end do
    end do
    !$omp end do nowait
  else
    !$omp single
    do j = iymin, iymax
      call interp_boundary(ctr(ixmin, j), ctr(ixmin, j), ctr(ixmin + 1, j), IDIRC)
      call interp_boundary(ctr(ixmax - 1, j), ctr(ixmax, j), ctr(ixmax, j), IDIRC)
    end do
    !$omp end single nowait
    !$omp do
    do j = iymin, iymax
      do i = ixmin + 1, ixmax - 1
        call interp_inner(ctr(i - 1, j), ctr(i, j), ctr(i + 1, j), IDIRC)
      end do
    end do
    !$omp end do nowait
  end if

  ! j direction
  !$omp single
  do i = ixmin, ixmax
    call interp_boundary(ctr(i, iymin), ctr(i, iymin), ctr(i, iymin + 1), JDIRC)
    call interp_boundary(ctr(i, iymax - 1), ctr(i, iymax), ctr(i, iymax), JDIRC)
  end do
  !$omp end single nowait

  !$omp do
  do j = iymin + 1, iymax - 1
    do i = ixmin, ixmax
      call interp_inner(ctr(i, j - 1), ctr(i, j), ctr(i, j + 1), JDIRC)
    end do
  end do
  !$omp end do nowait

  !$omp end parallel

end subroutine interpolation
! ------------------------------------------------------------
	! one-sided interpolation of the boundary cell
	subroutine interp_boundary(cell_L,cell_N,cell_R,DIRC)
    type(cell_center),intent(inout):: cell_N
    type(cell_center),intent(inout):: cell_L,cell_R
    integer,intent(in):: DIRC

    cell_N%sh(:,:,DIRC)=(cell_R%h-cell_L%h)/(0.5*cell_R%length(DIRC)+0.5*cell_L%length(DIRC))
    cell_N%sb(:,:,DIRC)=(cell_R%b-cell_L%b)/(0.5*cell_R%length(DIRC)+0.5*cell_L%length(DIRC))
    
    !Conservative variable slope
    cell_N%sw(:,DIRC)=(cell_R%w-cell_L%w)/(0.5*cell_R%length(DIRC)+0.5*cell_L%length(DIRC))
	end subroutine interp_boundary
! ------------------------------------------------------------
	! interpolation of the inner cell
	! using Van-Leer limiter
	subroutine interp_inner(cell_L,cell_N,cell_R,DIRC)
    type(cell_center),intent(in):: cell_L,cell_R
    type(cell_center),intent(inout):: cell_N
    integer,intent(in):: DIRC
    real(kind=8):: sL(unum,vnum),sR(unum,vnum)
    real(kind=8):: swL(DIM+2),swR(DIM+2)
    
    sL=(cell_N%h-cell_L%h)/(0.5*cell_N%length(DIRC)+0.5*cell_L%length(DIRC))
    sR=(cell_R%h-cell_N%h)/(0.5*cell_R%length(DIRC)+0.5*cell_N%length(DIRC))
    cell_N%sh(:,:,DIRC)=(sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
    
    sL=(cell_N%b-cell_L%b)/(0.5*cell_N%length(DIRC)+0.5*cell_L%length(DIRC))
    sR=(cell_R%b-cell_N%b)/(0.5*cell_R%length(DIRC)+0.5*cell_N%length(DIRC))
    cell_N%sb(:,:,DIRC)=(sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
    
    !Conservative variable slope
    swL=(cell_N%w-cell_L%w)/(0.5*cell_N%length(DIRC)+0.5*cell_L%length(DIRC))
    swR=(cell_R%w-cell_N%w)/(0.5*cell_R%length(DIRC)+0.5*cell_N%length(DIRC))
    cell_N%sw(:,DIRC)=(sign(UP,swR)+sign(UP,swL))*abs(swR)*abs(swL)/(abs(swR)+abs(swL)+smv)

	end subroutine interp_inner
! ------------------------------------------------------------
	subroutine evolution()
  use omp_lib
  integer :: i, j
  real(kind=8) :: bcwall(DIM+2)

  ! Set the number of threads for this parallel region
  call omp_set_num_threads(8)

  select case(case_id)
  !----------------------------------------------------------
  case(6)  ! Rayleigh-Taylor quarter-domain
  !   left/bottom: symmetry axes; right/top: impermeable specular walls.
  !----------------------------------------------------------
  !$omp parallel default(none) private(i, j) shared(ixmin, ixmax, iymin, iymax, vface, hface, ctr)
  !$omp single
  do j = iymin, iymax
    call calc_flux_specular_vertical(vface(ixmin, j), ctr(ixmin, j), IDIRC, RN)
    call calc_flux_specular_vertical(vface(ixmax+1, j), ctr(ixmax, j), IDIRC, RY)
  end do
  !$omp end single nowait
  !$omp do
  do j = iymin, iymax
    do i = ixmin + 1, ixmax
      call calc_flux_mtforcev(ctr(i - 1, j), vface(i, j), ctr(i, j), IDIRC)
    end do
  end do
  !$omp end do nowait
  !$omp single
  do i = ixmin, ixmax
    call calc_flux_specular_horizon(hface(i, iymin), ctr(i, iymin), JDIRC, RN)
    call calc_flux_specular_horizon(hface(i, iymax+1), ctr(i, iymax), JDIRC, RY)
  end do
  !$omp end single nowait
  !$omp do
  do j = iymin + 1, iymax
    do i = ixmin, ixmax
      call calc_flux_mtforceh(ctr(i, j - 1), hface(i, j), ctr(i, j), JDIRC)
    end do
  end do
  !$omp end do nowait
  !$omp end parallel

  !----------------------------------------------------------
  case(5)  ! Poiseuille: x-periodic, y-diffuse walls
  !   x direction: periodic ghost-cell copy, then internal flux loop
  !   spans ixmin..ixmax+1 so vface(ixmax+1) gets a flux too.
  !   y direction: diffuse isothermal wall (calc_flux_boundary) at
  !   iymin and iymax+1; prim=(rho,U,V,lambda), rho set by -R1/R2.
  !----------------------------------------------------------
  call bc_periodic_2d(IDIRC)
  !$omp parallel default(none) private(i, j, bcwall) shared(ixmin, ixmax, iymin, iymax, vface, hface, ctr, TwL, TwR)
  ! x direction: internal fluxes including the periodic boundary face
  !$omp do
  do j = iymin, iymax
    do i = ixmin, ixmax + 1
      call calc_flux_mtforcev(ctr(i - 1, j), vface(i, j), ctr(i, j), IDIRC)
    end do
  end do
  !$omp end do nowait
  ! y direction: diffuse isothermal walls
  !$omp do
  do i = ixmin, ixmax
    bcwall = (/0.d0, 0.d0, 0.d0, 1.d0/TwL/)
    call calc_flux_boundary(bcwall, hface(i, iymin), ctr(i, iymin), JDIRC, RN)
    bcwall = (/0.d0, 0.d0, 0.d0, 1.d0/TwR/)
    call calc_flux_boundary(bcwall, hface(i, iymax+1), ctr(i, iymax), JDIRC, RY)
  end do
  !$omp end do nowait
  !$omp do
  do j = iymin + 1, iymax
    do i = ixmin, ixmax
      call calc_flux_mtforceh(ctr(i, j - 1), hface(i, j), ctr(i, j), JDIRC)
    end do
  end do
  !$omp end do nowait
  !$omp end parallel

  case default
    write(*,*) "evolution: unknown case_id ",case_id
    stop
  end select

end subroutine evolution
! ------------------------------------------------------------
	subroutine update()
	!t=t^n
    real(kind=8),allocatable,dimension(:,:):: H_old,B_old !equilibrium distribution at t=t^n
    real(kind=8):: w_old(DIM+2)
    real(kind=8):: prim_old(DIM+2)
    real(kind=8):: tau_old
    
    !t=t^n+1
    real(kind=8),allocatable,dimension(:,:):: H,B !equilibrium distribution at t=t^{n+1}
    real(kind=8),allocatable,dimension(:,:):: H_plus,B_plus !Shakhov part
    
    real(kind=8):: prim(DIM+2) !primary variables at t^n and t^{n+1}
    real(kind=8):: qf(DIM)
    real(kind=8):: tau !collision time and t^n and t^{n+1}
    real(kind=8):: sum_res(DIM+2),mass_scale
    
    integer:: i,j,order
    real(kind=8)::error,errort
    
    real(kind=8),allocatable,dimension(:,:):: hsource
    real(kind=8),allocatable,dimension(:,:):: bsource
    REAL*8, ALLOCATABLE,dimension(:,:) :: ft1
    REAL*8, ALLOCATABLE ,dimension(:,:):: ft2
    
    real(kind=8):: fx,fy
    ALLOCATE(ft1(unum,vnum))
    ALLOCATE(ft2(unum,vnum))
    !allocate arrays
    allocate(H_old(unum,vnum))
    allocate(B_old(unum,vnum))
    allocate(H(unum,vnum))
    allocate(B(unum,vnum))
    allocate(H_plus(unum,vnum))
    allocate(B_plus(unum,vnum))
    allocate(hsource(unum,vnum))
    allocate(bsource(unum,vnum))

    !set initial value
    res=0.0
    sum_res=0.0
    mass_scale=0.0

    ! All physical cells are updated.  Boundary conditions act through
    ! boundary fluxes or true ghost cells.
    do j=iymin,iymax
      do i=ixmin,ixmax
            if (case_id==6) then
              fx=phi*cos(ctr(i,j)%theta)
              fy=phi*sin(ctr(i,j)%theta)
            else
              fx=phi
              fy=0.d0
            end if

            !--------------------------------------------------
            ! store W^n and calculate H^n,B^n,\tau^n
            !--------------------------------------------------
            w_old=ctr(i,j)%w !store W^n
            
            prim_old=get_primary(w_old) !convert to primary variables
            call reduced_maxwell(H_old,B_old,uspace,vspace,prim_old) !calculate Maxwellian
            tau_old=get_tau(prim_old) !calculate collision time \tau^n

            !--------------------------------------------------
            ! update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
            !--------------------------------------------------
            ctr(i,j)%w=ctr(i,j)%w+(vface(i,j)%flux-vface(i+1,j)%flux+hface(i,j)%flux-hface(i,j+1)%flux)/ctr(i,j)%area

            prim=get_primary(ctr(i,j)%w)

            !projection
            prim(2)=prim(2)+fx*dt
            prim(3)=prim(3)+fy*dt
            ctr(i,j)%w=get_conserved(prim)

            call reduced_maxwell(H,B,uspace,vspace,prim)
            tau=get_tau(prim)

            !--------------------------------------------------
            ! record residual
            !--------------------------------------------------
            sum_res=sum_res+(w_old-ctr(i,j)%w)**2
            ! Use one non-vanishing reference scale for every conserved
            ! component.  Component-wise normalization is singular for the
            ! Poiseuille wall-normal momentum, whose correct steady value is
            ! zero, and can therefore prevent a converged case from stopping.
            mass_scale=mass_scale+abs(ctr(i,j)%w(1))

            !--------------------------------------------------
            ! Shakhov part
            !--------------------------------------------------
            !heat flux at t=t^n
            qf=get_heat_flux(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim_old) 

            !h^+ = H+H^+ at t=t^n
            call shakhov_part(H_old,B_old,uspace,vspace,qf,prim_old,H_plus,B_plus) !H^+ and B^+
            H_old=H_old+H_plus !h^+
            B_old=B_old+B_plus !b^+

            !h^+ = H+H^+ at t=t^{n+1}
            call shakhov_part(H,B,uspace,vspace,qf,prim,H_plus,B_plus)
            H=H+H_plus
            B=B+B_plus

            !--------------------------------------------------
            ! body-force source term: finite-difference or Hermite spectral
            !--------------------------------------------------
            select case(trim(force_type))
            case('finite_difference')
                call cal_fsource1(ctr(i,j),vface(i,j))
                hsource=duh*dt*fx+dvh*dt*fy
                bsource=dub*dt*fx+dvb*dt*fy
            case('hermite')
                order=Nm
                call new_expansion_force_term(ctr(i,j)%h,ctr(i,j)%b, order, &
                    prim_old, fx, fy, ft1, ft2, error, errort)
                hsource=ft1*dt
                bsource=ft2*dt
            case('hermite_symmetric')
                order=Nm
                call symmetric_hermite_force_term_2d(ctr(i,j)%h,ctr(i,j)%b, &
                    order,prim_old,fx,fy,ft1,ft2,error,errort)
                hsource=ft1*dt
                bsource=ft2*dt
            case('hermite_symmetric_b')
                order=Nm
                call symmetric_b_hermite_force_term_2d(ctr(i,j)%h,ctr(i,j)%b, &
                    order,prim_old,fx,fy,ft1,ft2,error,errort)
                hsource=ft1*dt
                bsource=ft2*dt
            case default
                write(*,*) "unknown force_type '",trim(force_type),"'"
                stop
            end select

            !--------------------------------------------------
            ! update distribution function
            !--------------------------------------------------
            ctr(i,j)%h=(ctr(i,j)%h+(vface(i,j)%flux_h-vface(i+1,j)%flux_h+&
                        hface(i,j)%flux_h-hface(i,j+1)%flux_h)/ctr(i,j)%area+&
                        0.5*dt*(H/tau+(H_old-ctr(i,j)%h)/tau_old)-hsource)/(1.0+0.5*dt/tau)
            ctr(i,j)%b=(ctr(i,j)%b+(vface(i,j)%flux_b-vface(i+1,j)%flux_b+&
                        hface(i,j)%flux_b-hface(i,j+1)%flux_b)/ctr(i,j)%area+&
                        0.5*dt*(B/tau+(B_old-ctr(i,j)%b)/tau_old)-bsource)/(1.0+0.5*dt/tau)
        end do
    end do
            
    !final residual
    res=sqrt(grid%numElem*sum_res)/(mass_scale+SMV)
   
	end subroutine update
! ------------------------------------------------------------

end module SOLVER
