module GASTHEORY
use ALLDATA
implicit none
contains
! ------------------------------------------------------------
    pure function get_lambda_from_temperature(temperature)
    real(kind=8),intent(in):: temperature
    real(kind=8):: get_lambda_from_temperature

    get_lambda_from_temperature=1.d0/temperature
    end function get_lambda_from_temperature
! ------------------------------------------------------------
    pure function get_temperature_from_lambda(lambda)
    real(kind=8),intent(in):: lambda
    real(kind=8):: get_temperature_from_lambda

    get_temperature_from_lambda=1.d0/lambda
    end function get_temperature_from_lambda
! ------------------------------------------------------------
    subroutine reduced_maxwell(H,B,prim)
    real(kind=8),dimension(:),intent(out):: H,B
    real(kind=8),intent(in):: prim(DIM+2)

    H = prim(1)*(prim(3)/PI)**(1.0/2.0)*exp(-prim(3)*(uspace-prim(2))**2)
    B = h*innerk/(2.0*prim(3))
    !print *, "H at index", prim(3)
    end subroutine reduced_maxwell
! ------------------------------------------------------------
    subroutine shakhov_part(H,B,qf,prim,H_plus,B_plus)
    real(kind=8),dimension(:),intent(in) :: H,B
    real(kind=8),intent(in) :: qf
    real(kind=8),intent(in) :: prim(DIM+2)
    real(kind=8),dimension(:),intent(out) :: H_plus,B_plus

    H_plus = 0.8*(1-prandtl)*prim(3)**2/prim(1)*&
             (uspace-prim(2))*qf*(2*prim(3)*(uspace-prim(2))**2+innerk-5)*H
    B_plus = 0.8*(1-prandtl)*prim(3)**2/prim(1)*&
             (uspace-prim(2))*qf*(2*prim(3)*(uspace-prim(2))**2+innerk-3)*B
    end subroutine shakhov_part
! ------------------------------------------------------------
    function get_conserved(prim)
    real(kind=8),intent(in) :: prim(DIM+2)
    real(kind=8) :: get_conserved(DIM+2)

    get_conserved(1) = prim(1)
    get_conserved(2) = prim(1)*prim(2)
    get_conserved(3) = 0.5d0*prim(1)/prim(3)/(gamma-1.0)+0.5*prim(1)*prim(2)**2

    end function get_conserved
! ------------------------------------------------------------       
    function get_primary(w)
    real(kind=8),intent(in):: w(DIM+2)
    real(kind=8):: get_primary(DIM+2)

    get_primary(1) = w(1)
    get_primary(2) = w(2)/w(1)
    get_primary(3) = 0.5*w(1)/(gamma-1.0)/(w(3)-0.5*w(2)**2/w(1))
    end function get_primary
! ------------------------------------------------------------ 
    function get_gamma(K)
    real(kind=8),intent(in):: K
    real(kind=8):: get_gamma

    get_gamma=(K+3.0d0)/(K+1.0d0)
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
    ! calculate collision time: VHS model
    function get_tau(prim)
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8):: get_tau
    
    get_tau=muref*2*prim(3)**(1-omega)/prim(1)
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
    !get_mu=0.3125d0*sqrt(Rc*Tref/PI)*rhoinf*mfp*sqrt(2.d0)*PI
    return
    end function get_mu
! ------------------------------------------------------------
    function get_heat_flux(h,b,prim)
    real(kind=8),dimension(:),intent(in):: h,b
    real(kind=8),intent(in):: prim(DIM+2)
    real(kind=8):: get_heat_flux
    
    get_heat_flux = 0.5d0*(sum(weight*(uspace-prim(2))*(uspace-prim(2))**2*h)+sum(weight*(uspace-prim(2))*b)) 
    end function get_heat_flux
! ------------------------------------------------------------
end module GASTHEORY
