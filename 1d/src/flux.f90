module FLUX
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
implicit none

integer,parameter :: MNUM = 6 !number of normal velocity moments
integer,parameter :: MTUM = 4 !number of tangential velocity moments

contains
! ------------------------------------------------------------
    subroutine calc_flux(cell_L,face,cell_R)
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell_L,cell_R

    !Heaviside step function
    integer,allocatable,dimension(:):: delta

    !interface variable
    real(kind=8),allocatable,dimension(:):: h,b
    real(kind=8),allocatable,dimension(:):: H0,B0
    real(kind=8),allocatable,dimension(:):: H_plus,B_plus
    real(kind=8),allocatable,dimension(:):: sh,sb 
    real(kind=8):: w(DIM+2),prim(DIM+2)
    real(kind=8):: qf
    real(kind=8):: sw(DIM+2)
    real(kind=8):: aL(DIM+2),aR(DIM+2),aT(DIM+2)

    !moments variable
    real(kind=8):: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
    real(kind=8):: Mau_0(3),Mau_L(3),Mau_R(3),Mau_T(3) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<A*u*\psi>
    real(kind=8):: tau
    real(kind=8):: Mt(6)

    !original velocity space
    real(kind=8),allocatable,dimension(:):: vg
    real(kind=8),allocatable,dimension(:):: vf

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    allocate(delta(unum))
    allocate(h(unum))
    allocate(b(unum))
    allocate(sh(unum))
    allocate(sb(unum))
    allocate(H0(unum))
    allocate(B0(unum))
    allocate(H_plus(unum))
    allocate(B_plus(unum))
    allocate(vg(unum))
    allocate(vf(unum))

    !original velocity space
    vg=uspace-0.25d0*phi*dt
    vf=uspace-0.5d0*phi*dt

    !Heaviside step function
    delta = (sign(UP,vf)+1.d0)/2.d0

    !--------------------------------------------------
    ! UPWIND reconstruction
    !--------------------------------------------------
    h = (cell_L%h+0.5*cell_L%length*cell_L%sh)*delta+&
        (cell_R%h-0.5*cell_R%length*cell_R%sh)*(1-delta)
    b = (cell_L%b+0.5*cell_L%length*cell_L%sb)*delta+&
        (cell_R%b-0.5*cell_R%length*cell_R%sb)*(1-delta)
    sh = cell_L%sh*delta+cell_R%sh*(1-delta)
    sb = cell_L%sb*delta+cell_R%sb*(1-delta)

    !--------------------------------------------------
    ! obtain macroscopic variables at interface
    !--------------------------------------------------
    !conservative variables W_0 
    w(1) = sum(weight*h)
    w(2) = sum(weight*uspace*h)
    w(3) = 0.5*(sum(weight*uspace**2*h)+sum(weight*b))

    !convert to primary variables
    prim = get_primary(w)

    !heat flux
    qf = get_heat_flux(h,b,prim) 

    !--------------------------------------------------
    ! calculate a^L,a^R
    !--------------------------------------------------
    sw = (w-cell_L%w)/(0.5*cell_L%length)
    aL = micro_slope(prim,sw)

    sw = (cell_R%w-w)/(0.5*cell_R%length)
    aR = micro_slope(prim,sw)

    !--------------------------------------------------
    ! calculate time slope of W and A
    !--------------------------------------------------
    !<u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
    call calc_moment_u(prim,Mu,Mxi,Mu_L,Mu_R) 

    Mau_L = moment_au(aL,Mu_L,Mxi,1) !<aL*u*\psi>_{>0}
    Mau_R = moment_au(aR,Mu_R,Mxi,1) !<aR*u*\psi>_{<0}

    sw = -prim(1)*(Mau_L+Mau_R) !time slope of W
    aT = micro_slope(prim,sw) !calculate A

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau=get_tau(prim)

    Mt(4)=tau*(1.0-exp(-dt/tau))
    Mt(5)=-tau*dt*exp(-dt/tau)+tau*Mt(4)
    Mt(1)=dt-Mt(4)
    Mt(2)=-tau*Mt(1)+Mt(5) 
    Mt(3)=dt**2/2.0-tau*Mt(1)
    Mt(6)=tau**3-0.5d0*exp(-dt/tau)*tau*(dt**2+2*dt*tau+2*tau**2)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Mau_0 = moment_uv(Mu,Mxi,1,0) !<u*\psi>
    Mau_L = moment_au(aL,Mu_L,Mxi,2) !<aL*u^2*\psi>_{>0}
    Mau_R = moment_au(aR,Mu_R,Mxi,2) !<aR*u^2*\psi>_{<0}
    Mau_T = moment_au(aT,Mu,Mxi,1) !<A*u*\psi>

    face%flux = Mt(1)*prim(1)*Mau_0+Mt(2)*prim(1)*(Mau_L+Mau_R)+Mt(3)*prim(1)*Mau_T

    !--------------------------------------------------
    ! gravity flux related to g0
    !--------------------------------------------------
    Mau_0=moment_uv(Mu,Mxi,1,0)
    Mau_L=moment_au(aL,Mu_L,Mxi,1)
    Mau_R=moment_au(aR,Mu_R,Mxi,1)

    face%flux=face%flux-tau*Mt(2)*prim(1)*phi*(Mau_L+Mau_R)-&
              Mt(6)*prim(1)*phi*(Mau_L+Mau_R)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    !Maxwellian distribution H0 and B0
    call reduced_maxwell(H0,B0,prim)

    !Shakhov part H+ and B+
    call shakhov_part(H0,B0,qf,prim,H_plus,B_plus)

    !macro flux related to g+ and f0
    face%flux(1) = face%flux(1)+Mt(1)*sum(weight*uspace*H_plus)+Mt(4)*sum(weight*uspace*h)-Mt(5)*sum(weight*uspace**2*sh)
    face%flux(2) = face%flux(2)+Mt(1)*sum(weight*uspace**2*H_plus)+Mt(4)*sum(weight*uspace**2*h)-Mt(5)*sum(weight*uspace**3*sh)
    face%flux(3) = face%flux(3)+&
                    Mt(1)*0.5*(sum(weight*uspace*uspace**2*H_plus)+sum(weight*uspace*B_plus))+&
                    Mt(4)*0.5*(sum(weight*uspace*uspace**2*h)+sum(weight*uspace*b))-&
                    Mt(5)*0.5*(sum(weight*uspace**2*uspace**2*sh)+sum(weight*uspace**2*sb))

    !gravity concern
    face%flux(1)=face%flux(1)+phi*Mt(6)*sum(weight*uspace*sh)
    face%flux(2)=face%flux(2)+phi*Mt(6)*sum(weight*uspace**2*sh)
    face%flux(3)=face%flux(3)+phi*Mt(6)*0.5*(sum(weight*(uspace**2)*uspace*sh)+sum(weight*uspace*sb))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    face%flux_h=Mt(1)*uspace*(H0+H_plus)+&
                Mt(2)*uspace**2*(aL(1)*H0+aL(2)*uspace*H0+0.5*aL(3)*(uspace**2*H0+B0))*delta+&
                Mt(2)*uspace**2*(aR(1)*H0+aR(2)*uspace*H0+0.5*aR(3)*(uspace**2*H0+B0))*(1-delta)+&
                Mt(3)*uspace*(aT(1)*H0+aT(2)*uspace*H0+0.5*aT(3)*(uspace**2*H0+B0))+&
                Mt(4)*uspace*h-Mt(5)*uspace**2*sh

    face%flux_b = Mt(1)*uspace*(B0+B_plus)+&
                  Mt(2)*uspace**2*(aL(1)*B0+aL(2)*uspace*B0+0.5*aL(3)*(uspace**2*B0+Mxi(2)*H0))*delta+&
                  Mt(2)*uspace**2*(aR(1)*B0+aR(2)*uspace*B0+0.5*aR(3)*(uspace**2*B0+Mxi(2)*H0))*(1-delta)+&
                  Mt(3)*uspace*(aT(1)*B0+aT(2)*uspace*B0+0.5*aT(3)*(uspace**2*B0+Mxi(2)*H0))+&
                  Mt(4)*uspace*b-Mt(5)*uspace**2*sb

    !gravity concern
    face%flux_h=face%flux_h-&
                tau*Mt(2)*phi*uspace*(aL(1)*H0+aL(2)*uspace*H0+0.5*aL(3)*((uspace**2)*H0+B0))*delta-&
                tau*Mt(2)*phi*uspace*(aR(1)*H0+aR(2)*uspace*H0+0.5*aR(3)*((uspace**2)*H0+B0))*(1-delta)+&
                Mt(6)*phi*uspace*sh
    face%flux_b=face%flux_b-&
                tau*Mt(2)*phi*uspace*(aL(1)*B0+aL(2)*uspace*B0+0.5*aL(3)*((uspace**2)*B0+Mxi(2)*H0))*delta-&
                tau*Mt(2)*phi*uspace*(aR(1)*B0+aR(2)*uspace*B0+0.5*aR(3)*((uspace**2)*B0+Mxi(2)*H0))*(1-delta)+&
                Mt(6)*phi*uspace*sb

    end subroutine calc_flux
! ------------------------------------------------------------
    subroutine calc_flux_specular(face,cell,rot) 
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell
    integer,intent(in):: rot
    
    !Heaviside step function
    integer,allocatable,dimension(:):: delta 
    
    !interface variable
    real(kind=8),allocatable,dimension(:):: h,b
    real(kind=8),allocatable,dimension(:):: H0,B0
    real(kind=8):: prim(DIM+2) !boundary condition in local frame
    real(kind=8),allocatable,dimension(:):: temph,tempb

    integer:: i

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !allocate array
    allocate(delta(unum))
    allocate(h(unum))
    allocate(b(unum))
    allocate(H0(unum))
    allocate(B0(unum))
    allocate(temph(unum))
    allocate(tempb(unum))
    
    delta=(sign(UP,uspace)*rot+1)/2

    h=cell%h-rot*0.5*cell%length*cell%sh
    b=cell%b-rot*0.5*cell%length*cell%sb

    do i=1,unum
        temph(i)=h(unum-i+1)
        tempb(i)=b(unum-i+1)
    end do
    h=temph*delta+h*(1-delta)
    b=tempb*delta+b*(1-delta)

    face%flux(1)=sum(weight*uspace*h)
    face%flux(2)=sum(weight*uspace*uspace*h)
    face%flux(3)=0.5*sum(weight*uspace*((uspace**2)*h+b))

    face%flux_h=uspace*h
    face%flux_b=uspace*b

    face%flux=dt*face%flux
    face%flux_h=dt*face%flux_h
    face%flux_b=dt*face%flux_b

    end subroutine calc_flux_specular
! ------------------------------------------------------------
subroutine calc_flux_adiabatic(face,cell,rot) 
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell
    integer,intent(in):: rot
    
    !Heaviside step function
    integer,allocatable,dimension(:):: delta 
    
    !interface variable
    real(kind=8),allocatable,dimension(:):: h,b
    real(kind=8),allocatable,dimension(:):: H0,B0
    real(kind=8):: prim(DIM+2) !boundary condition in local frame
    real(kind=8):: R1,R2
    real(kind=8):: T1,T2

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !allocate array
    allocate(delta(unum))
    allocate(h(unum))
    allocate(b(unum))
    allocate(H0(unum))
    allocate(B0(unum))
    
    delta=(sign(UP,uspace)*rot+1)/2

    h=cell%h-rot*0.5*cell%length*cell%sh
    b=cell%b-rot*0.5*cell%length*cell%sb

    T1=sum(weight*(uspace**3*h+uspace*b)*(1-delta))
    T2=sum(weight*uspace*h*(1-delta))
    prim(3)=T2/T1
    write(*,*) "lambda:",prim(3)
    !R1=sum(weight*uspace*h*(1-delta))
    !R2=(prim(3)/PI)*sum(weight*uspace*exp(-prim(4)*((uspace-prim(2))**2))*delta)
    !prim(1)=-R1/R2
    prim(1)=2.d0*dsqrt(PI)*dsqrt(prim(3))*sum(weight*uspace*h*(1-delta))

    prim(2)=0
    call reduced_maxwell(H0,B0,prim)

    h=H0*delta+h*(1-delta)
    b=B0*delta+b*(1-delta)

    face%flux(1)=sum(weight*uspace*h)
    face%flux(2)=sum(weight*uspace*uspace*h)
    face%flux(3)=0.5*sum(weight*uspace*((uspace**2)*h+b))

    face%flux_h=uspace*h
    face%flux_b=uspace*b

    face%flux=dt*face%flux
    face%flux_h=dt*face%flux_h
    face%flux_b=dt*face%flux_b

    end subroutine calc_flux_adiabatic
! ------------------------------------------------------------ 
    subroutine calc_flux_boundary(bc,face,cell,rot) 
    real(kind=8),intent(in):: bc(DIM+2)
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell
    integer,intent(in):: rot
    
    !Heaviside step function
    integer,allocatable,dimension(:):: delta 
    
    !interface variable
    real(kind=8),allocatable,dimension(:):: h,b
    real(kind=8),allocatable,dimension(:):: H0,B0
    real(kind=8):: prim(DIM+2) !boundary condition in local frame
    real(kind=8):: R1,R2

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !allocate array
    allocate(delta(unum))
    allocate(h(unum))
    allocate(b(unum))
    allocate(H0(unum))
    allocate(B0(unum))
    
    delta=(sign(UP,uspace)*rot+1)/2

    prim=bc

    h=cell%h-rot*0.5*cell%length*cell%sh
    b=cell%b-rot*0.5*cell%length*cell%sb

    R1=sum(weight*uspace*h*(1-delta))
    ! Construct the unit-density wall Maxwellian with the same routine used
    ! everywhere else.  Its discrete outgoing half-space moment is the exact
    ! coefficient needed to balance the incoming discrete mass flux.
    prim(1)=1.d0
    call reduced_maxwell(H0,B0,prim)
    R2=sum(weight*uspace*H0*delta)
    if (abs(R2) <= tiny(1.d0)) then
      write(*,*) 'calc_flux_boundary: zero outgoing wall mass coefficient'
      error stop
    end if
    prim(1)=-R1/R2

    call reduced_maxwell(H0,B0,prim)

    h=H0*delta+h*(1-delta)
    b=B0*delta+b*(1-delta)

    face%flux(1)=sum(weight*uspace*h)
    face%flux(2)=sum(weight*uspace*uspace*h)
    face%flux(3)=0.5*sum(weight*uspace*((uspace**2)*h+b))

    face%flux_h=uspace*h
    face%flux_b=uspace*b

    face%flux=dt*face%flux
    face%flux_h=dt*face%flux_h
    face%flux_b=dt*face%flux_b

    end subroutine calc_flux_boundary
! ------------------------------------------------------------
    function micro_slope(prim,sw)
        real(kind=8),intent(in) :: prim(3),sw(3)
        real(kind=8) :: micro_slope(3)

        micro_slope(3) = 4.0*prim(3)**2/(innerk+1)/prim(1)*(2.0*sw(3)-2.0*prim(2)*sw(2)+sw(1)*(prim(2)**2-0.5*(innerk+1)/prim(3)))

        micro_slope(2) = 2.0*prim(3)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*micro_slope(3)
        micro_slope(1) = sw(1)/prim(1)-prim(2)*micro_slope(2)-0.5*(prim(2)**2+0.5*(innerk+1)/prim(3))*micro_slope(3)
    end function micro_slope
! ------------------------------------------------------------
    subroutine calc_moment_u(prim,Mu,Mxi,Mu_L,Mu_R)
        real(kind=8),intent(in) :: prim(3)
        real(kind=8),intent(out) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)
        real(kind=8),intent(out) :: Mxi(0:2)
        integer :: i

        !moments of normal velocity
        Mu_L(0) = 0.5*erfc(-sqrt(prim(3))*prim(2))
        Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))
        Mu_R(0) = 0.5*erfc(sqrt(prim(3))*prim(2))
        Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))

        do i=2,MNUM
            Mu_L(i) = prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)/prim(3)
            Mu_R(i) = prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)/prim(3)
        end do

        Mu = Mu_L+Mu_R

        !moments of \xi
        Mxi(0) = 1.0 !<\xi^0>
        Mxi(1) = 0.5*innerk/prim(3) !<\xi^2>
        Mxi(2) = (innerk**2+2.0*innerk)/(4.0*prim(3)**2) !<\xi^4>
    end subroutine calc_moment_u
! ------------------------------------------------------------
    function moment_uv(Mu,Mxi,alpha,delta)
        real(kind=8),intent(in) :: Mu(0:MNUM),Mxi(0:2)
        integer,intent(in) :: alpha,delta
        real(kind=8) :: moment_uv(3)

        moment_uv(1) = Mu(alpha)*Mxi(delta/2)
        moment_uv(2) = Mu(alpha+1)*Mxi(delta/2)
        moment_uv(3) = 0.5*(Mu(alpha+2)*Mxi(delta/2)+Mu(alpha)*Mxi((delta+2)/2))
    end function moment_uv
! ------------------------------------------------------------
    function moment_au(a,Mu,Mxi,alpha)
        real(kind=8),intent(in) :: a(3)
        real(kind=8),intent(in) :: Mu(0:MNUM),Mxi(0:2)
        integer,intent(in) :: alpha
        real(kind=8) :: moment_au(3)

        moment_au = a(1)*moment_uv(Mu,Mxi,alpha+0,0)+&
                    a(2)*moment_uv(Mu,Mxi,alpha+1,0)+&
                    0.5*a(3)*moment_uv(Mu,Mxi,alpha+2,0)+&
                    0.5*a(3)*moment_uv(Mu,Mxi,alpha+0,2)
    end function moment_au
! ------------------------------------------------------------
end module FLUX
