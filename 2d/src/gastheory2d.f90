module GASTHEORY
use ALLDATA
implicit none
contains
! ------------------------------------------------------------
    subroutine reduced_maxwell(H,B,vn,vt,prim)
    real(kind=8),dimension(:,:),intent(out):: H,B
    real(kind=8),dimension(:,:),intent(in):: vn,vt
    real(kind=8),dimension(:,:),intent(in):: prim(DIM+2)

    H=prim(1)*(prim(4)/PI)*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))
    B=H*innerk/(2.0*prim(4))
    end subroutine reduced_maxwell
! ------------------------------------------------------------
    subroutine shakhov_part(H,B,vn,vt,qf,prim,H_plus,B_plus)
    real(kind=8),dimension(:,:),intent(in):: H,B
    real(kind=8),dimension(:,:),intent(in):: vn,vt
    real(kind=8),intent(in):: qf(DIM)
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8),dimension(:,:),intent(out):: H_plus,B_plus

    H_plus=0.8*(1-prandtl)*prim(4)**2/prim(1)*&
            ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+innerk-5)*H
    B_plus=0.8*(1-prandtl)*prim(4)**2/prim(1)*&
            ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+innerk-3)*B
    end subroutine shakhov_part
! ------------------------------------------------------------
    function get_conserved(prim)
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8):: get_conserved(DIM+2)

    get_conserved(1)=prim(1)
    get_conserved(2)=prim(1)*prim(2)
    get_conserved(3)=prim(1)*prim(3)
    get_conserved(4)=0.5*prim(1)/prim(4)/(gamma-1.0)+0.5*prim(1)*(prim(2)**2+prim(3)**2)
    end function get_conserved
! ------------------------------------------------------------       
    function get_primary(w)
    real(kind=8),intent(in):: w(DIM+2)
    real(kind=8):: get_primary(DIM+2)

    get_primary(1)=w(1)
    get_primary(2)=w(2)/w(1)
    get_primary(3)=w(3)/w(1)
    get_primary(4)=0.5*w(1)/(gamma-1.0)/(w(4)-0.5*(w(2)**2+w(3)**2)/w(1))
    end function get_primary
! ------------------------------------------------------------   
    function global_frame(w,cosa,sina)
    real(kind=8),intent(in):: w(DIM+2)
    real(kind=8),intent(in):: cosa,sina
    real(kind=8):: global_frame(DIM+2)

    global_frame(1)=w(1)
    global_frame(2)=w(2)*cosa-w(3)*sina
    global_frame(3)=w(2)*sina+w(3)*cosa
    global_frame(4)=w(4)
    end function global_frame
! ------------------------------------------------------------ 
    function local_frame(w,cosa,sina)
    real(kind=8),intent(in):: w(DIM+2)
    real(kind=8),intent(in):: cosa,sina
    real(kind=8):: local_frame(DIM+2)

    local_frame(1)=w(1)
    local_frame(2)=w(2)*cosa+w(3)*sina  !u'=ucosa+vsina  
    local_frame(3)=w(3)*cosa-w(2)*sina  !v'=u(-sina)+vcosa
    local_frame(4)=w(4)
    end function local_frame
! ------------------------------------------------------------ 
    function get_gamma(K)
    real(kind=8),intent(in):: K
    real(kind=8):: get_gamma

    get_gamma=(K+4.0d0)/(K+2.0d0)
    return
    end function get_gamma
! ------------------------------------------------------------
    function get_sos(lam)
    real(kind=8),intent(in):: lam
    real(kind=8):: get_sos 

    get_sos=sqrt(0.5d0*gamma/lam)
    return
    end function get_sos
! ------------------------------------------------------------
    ! Calculate the collision time with the model used by each reference case.
    function get_tau(prim)
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8):: get_tau
    
    if (case_id==5) then
      get_tau=muref/prim(1)
    else
      get_tau=muref*2*prim(4)**(1-omega)/prim(1)
    end if
    return
    end function get_tau
! ------------------------------------------------------------
    ! calculate nondimensionalized viscosity coefficient
    ! alpha,omega :index related to HS/VHS/VSS model
    function get_mu(kn,alpha,omega)
    real(kind=8),intent(in):: kn,alpha,omega
    real(kind=8):: get_mu
        
    ! Hard Sphere model(HS model)
    get_mu=5.0d0*(alpha+1.0d0)*(alpha+2.0d0)*sqrt(PI)/&
    (4*alpha*(5-2*omega)*(7-2*omega))*kn !Kn  
    !get_mu=0.3125d0*dsqrt(Rc*Tref/PI)*rhoinf*mfp*dsqrt(2.d0)*PI

    return
    end function get_mu
! ------------------------------------------------------------
    function suther(temperature,mu)
    real(kind=8),intent(in):: temperature,mu
    real(kind=8):: suther

    suther = mu*(temperature/Tref)**omega

    end function suther
! ------------------------------------------------------------
    function get_heat_flux(h,b,vn,vt,prim)
    real(kind=8),dimension(:,:),intent(in):: h,b
    real(kind=8),dimension(:,:),intent(in):: vn,vt
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8):: get_heat_flux(DIM) !heat flux in normal and tangential direction

    get_heat_flux(1)=0.5*(sum(weight*(vn-prim(2))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vn-prim(2))*b)) 
    get_heat_flux(2)=0.5*(sum(weight*(vt-prim(3))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vt-prim(3))*b)) 
    end function get_heat_flux
! ------------------------------------------------------------
    function get_error(h,b,vn,vt,prim)
    real(kind=8),dimension(:,:),intent(in):: h,b
    real(kind=8),dimension(:,:),intent(in):: vn,vt
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8):: get_error(1) !heat flux in normal and tangential direction

    get_error(1)=1!(h-prim(1)*(prim(4)/PI)*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)))/h 
    end function get_error
    ! ------------------------------------------------------------
end module GASTHEORY
