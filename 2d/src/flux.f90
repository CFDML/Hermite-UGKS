module FLUX
use ALLDATA
use GRIDMAKER
use GASTHEORY
implicit none

integer,parameter :: MNUM = 6 !number of normal velocity moments
integer,parameter :: MTUM = 6 !number of tangential velocity moments

contains
! ------------------------------------------------------------
    subroutine face_acceleration(face,fn,ft)
    type(cell_interface),intent(in):: face
    real(kind=8),intent(out):: fn,ft
    real(kind=8):: ax,ay

    if (case_id==5) then
      ax=phi
      ay=0.d0
    else
      ax=phi*cos(face%theta)
      ay=phi*sin(face%theta)
    end if
    fn=ax*face%cosa+ay*face%sina
    ft=ay*face%cosa-ax*face%sina
    end subroutine face_acceleration
! ------------------------------------------------------------
    subroutine calc_flux_mtforcev(cell_L,face,cell_R,dirc)
    type(cell_interface),intent(inout):: face
    integer,intent(in):: dirc
    type(cell_center),intent(in):: cell_L,cell_R

    !Heaviside step function
    integer:: delta(unum,vnum)

    !Local micro velocity space
    real(kind=8):: vn(unum,vnum),vt(unum,vnum)

    !interface variable
    real(kind=8):: h(unum,vnum),b(unum,vnum)
    real(kind=8):: H0(unum,vnum),B0(unum,vnum)
    real(kind=8):: H_plus(unum,vnum),B_plus(unum,vnum)
    real(kind=8):: shn(unum,vnum),sbn(unum,vnum),sht(unum,vnum),sbt(unum,vnum)
    real(kind=8):: w(4),prim(4)
    real(kind=8):: qf(2)
    real(kind=8):: sw(4)
    real(kind=8):: aL(4),aR(4),aT(4),ab(4)

    !moments variable
    real(kind=8):: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
    real(kind=8):: Mau_0(4) !<u\psi>
    real(kind=8):: Mau_L(4),Mau_R(4) !<aL*u^n*\psi>,<aR*u^n*\psi>
    real(kind=8):: Mbu(4) !<b*u^n*\psi>
    real(kind=8):: Mau_T(4) !<A*u*\psi>
    real(kind=8):: tau
    real(kind=8):: Mt(6)

    !trace local velocity space
    real(kind=8):: vf(unum,vnum)
    real(kind=8):: fn,ft

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !--------------------------------------------------
    ! local variable
    !--------------------------------------------------
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina

    call face_acceleration(face,fn,ft)
    vf=vn-0.5d0*fn*dt

    delta=(sign(UP,vf)+1.d0)/2.d0

    !--------------------------------------------------
    ! reconstruct initial distribution
    !--------------------------------------------------
    ! UPWIND reconstruction
    h=(cell_L%h+0.5*cell_L%length(DIRC)*cell_L%sh(:,:,DIRC))*delta+&
      (cell_R%h-0.5*cell_R%length(DIRC)*cell_R%sh(:,:,DIRC))*(1-delta)
    b=(cell_L%b+0.5*cell_L%length(DIRC)*cell_L%sb(:,:,DIRC))*delta+&
      (cell_R%b-0.5*cell_R%length(DIRC)*cell_R%sb(:,:,DIRC))*(1-delta)
    shn=cell_L%sh(:,:,DIRC)*delta+cell_R%sh(:,:,DIRC)*(1-delta)
    sht=cell_L%sh(:,:,mod(DIRC,2)+1)*delta+cell_R%sh(:,:,mod(DIRC,2)+1)*(1-delta)
    sbn=cell_L%sb(:,:,DIRC)*delta+cell_R%sb(:,:,DIRC)*(1-delta)
    sbt=cell_L%sb(:,:,mod(DIRC,2)+1)*delta+cell_R%sb(:,:,mod(DIRC,2)+1)*(1-delta)

    !--------------------------------------------------
    ! obtain macroscopic variables at interface
    !--------------------------------------------------
    !conservative variables W_0 
    w(1)=sum(weight*h)
    w(2)=sum(weight*vn*h)
    w(3)=sum(weight*vt*h)
    w(4)=0.5*(sum(weight*(vn**2+vt**2)*h)+sum(weight*b))

    !convert to primary variables
    prim=get_primary(w)

    !heat flux
    qf=get_heat_flux(h,b,vn,vt,prim) 

    !--------------------------------------------------
    ! calculate a^L,a^R
    !--------------------------------------------------
    sw=(w-local_frame(cell_L%w,face%cosa,face%sina))/(0.5*cell_L%length(DIRC))
    aL=micro_slope(prim,sw)

    sw=(local_frame(cell_R%w,face%cosa,face%sina)-w)/(0.5*cell_R%length(DIRC))
    aR=micro_slope(prim,sw)

    !--------------------------------------------------
    ! calculate b
    !--------------------------------------------------
    sw(1)=sum(weight*sht)
    sw(2)=sum(weight*vn*sht)
    sw(3)=sum(weight*vt*sht)
    sw(4)=0.5d0*(sum(weight*(vn**2+vt**2)*sht)+sum(weight*sbt))
    ab=micro_slope(prim,sw)

    !--------------------------------------------------
    ! calculate time slope of W and A
    !--------------------------------------------------
    !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
    call calc_moment(prim,Mu,Mv,Mxi,Mu_L,Mu_R) 

    Mau_L=moment_au(aL,Mu_L,Mv,Mxi,1,0) !<aL*u*\psi>_{>0}
    Mau_R=moment_au(aR,Mu_R,Mv,Mxi,1,0) !<aR*u*\psi>_{<0}
    Mbu = moment_au(ab,Mu,Mv,Mxi,0,1)

    sw=-prim(1)*(Mau_L+Mau_R+Mbu) !time slope of W
    call reduced_maxwell(H0,B0,vn,vt,prim)
    call add_force_compatibility(sw,prim,vn,vt,H0,B0,fn,ft)
    aT=micro_slope(prim,sw) !calculate A

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
    Mau_0=moment_uv(Mu,Mv,Mxi,1,0,0) !<u*\psi>
    Mau_L=moment_au(aL,Mu_L,Mv,Mxi,2,0) !<aL*u^2*\psi>_{>0}
    Mau_R=moment_au(aR,Mu_R,Mv,Mxi,2,0) !<aR*u^2*\psi>_{<0}
    Mbu=moment_au(ab,Mu,Mv,Mxi,1,1) !!<b*u*v*\psi>
    Mau_T=moment_au(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>

    face%flux=Mt(1)*prim(1)*Mau_0+Mt(2)*prim(1)*(Mau_L+Mau_R+Mbu)+Mt(3)*prim(1)*Mau_T

    !--------------------------------------------------
    ! gravity flux related to g0
    !--------------------------------------------------
    !normal direction
    Mau_L=moment_au(aL,Mu_L,Mv,Mxi,1,0)
    Mau_R=moment_au(aR,Mu_R,Mv,Mxi,1,0)

    face%flux=face%flux-tau*Mt(2)*prim(1)*fn*(Mau_L+Mau_R)-&
              Mt(6)*prim(1)*fn*(Mau_L+Mau_R)
    !tangential direction
    Mbu=moment_au(ab,Mu,Mv,Mxi,0,1)
    face%flux=face%flux-tau*Mt(2)*prim(1)*ft*Mbu-&
              Mt(6)*prim(1)*ft*Mbu

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    !Maxwellian distribution H0 and B0
    call reduced_maxwell(H0,B0,vn,vt,prim)

    !Shakhov part H+ and B+
    call shakhov_part(H0,B0,vn,vt,qf,prim,H_plus,B_plus)

    !macro flux related to g+ and f0
    face%flux(1)=face%flux(1)+Mt(1)*sum(weight*vn*H_plus)+Mt(4)*sum(weight*vn*h)-Mt(5)*sum(weight*vn**2*shn)-Mt(5)*sum(weight*vt*vn*sht)
    face%flux(2)=face%flux(2)+Mt(1)*sum(weight*vn*vn*H_plus)+Mt(4)*sum(weight*vn*vn*h)-Mt(5)*sum(weight*vn*vn**2*shn)-Mt(5)*sum(weight*vn*vt*vn*sht)
    face%flux(3)=face%flux(3)+Mt(1)*sum(weight*vt*vn*H_plus)+Mt(4)*sum(weight*vt*vn*h)-Mt(5)*sum(weight*vt*vn**2*shn)-Mt(5)*sum(weight*vt*vt*vn*sht)
    face%flux(4)=face%flux(4)+&
                 Mt(1)*0.5*(sum(weight*vn*(vn**2+vt**2)*H_plus)+sum(weight*vn*B_plus))+&
                 Mt(4)*0.5*(sum(weight*vn*(vn**2+vt**2)*h)+sum(weight*vn*b))-&
                 Mt(5)*0.5*(sum(weight*vn**2*(vn**2+vt**2)*shn)+sum(weight*vn**2*sbn))-&
                 Mt(5)*0.5*(sum(weight*vt*vn*(vn**2+vt**2)*sht)+sum(weight*vt*vn*sbt))

    !gravity related to g+ and f0
    face%flux(1)=face%flux(1)+fn*Mt(6)*sum(weight*vn*shn)
    face%flux(2)=face%flux(2)+fn*Mt(6)*sum(weight*vn*vn*shn)
    face%flux(3)=face%flux(3)+fn*Mt(6)*sum(weight*vt*vn*shn)
    face%flux(4)=face%flux(4)+fn*Mt(6)*0.5*(sum(weight*(vn**2+vt**2)*vn*shn)+sum(weight*vn*sbn))

    face%flux(1)=face%flux(1)+ft*Mt(6)*sum(weight*vn*sht)
    face%flux(2)=face%flux(2)+ft*Mt(6)*sum(weight*vn*vn*sht)
    face%flux(3)=face%flux(3)+ft*Mt(6)*sum(weight*vt*vn*sht)
    face%flux(4)=face%flux(4)+ft*Mt(6)*0.5*(sum(weight*(vn**2+vt**2)*vn*sht)+sum(weight*vn*sbt))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    face%flux_h=Mt(1)*vn*(H0+H_plus)+&
                Mt(2)*vn**2*(aL(1)*H0+aL(2)*vn*H0+aL(3)*vt*H0+0.5*aL(4)*((vn**2+vt**2)*H0+B0))*delta+& 
                Mt(2)*vn**2*(aR(1)*H0+aR(2)*vn*H0+aR(3)*vt*H0+0.5*aR(4)*((vn**2+vt**2)*H0+B0))*(1-delta)+&
                Mt(2)*vn*vt*(ab(1)*H0+ab(2)*vn*H0+ab(3)*vt*H0+0.5*ab(4)*((vn**2+vt**2)*H0+B0))+&
                Mt(3)*vn*(aT(1)*H0+aT(2)*vn*H0+aT(3)*vt*H0+0.5*aT(4)*((vn**2+vt**2)*H0+B0))+&
                Mt(4)*vn*h-Mt(5)*vn**2*shn-Mt(5)*vn*vt*sht

    face%flux_b=Mt(1)*vn*(B0+B_plus)+&
                Mt(2)*vn**2*(aL(1)*B0+aL(2)*vn*B0+aL(3)*vt*B0+0.5*aL(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*delta+&
                Mt(2)*vn**2*(aR(1)*B0+aR(2)*vn*B0+aR(3)*vt*B0+0.5*aR(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*(1-delta)+&
                Mt(2)*vn*vt*(ab(1)*B0+ab(2)*vn*B0+ab(3)*vt*B0+0.5*ab(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                Mt(3)*vn*(aT(1)*B0+aT(2)*vn*B0+aT(3)*vt*B0+0.5*aT(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                Mt(4)*vn*b-Mt(5)*vn**2*sbn-Mt(5)*vn*vt*sbt

    !gravity for distribution function
    face%flux_h=face%flux_h-&
                tau*Mt(2)*fn*vn*(aL(1)*H0+aL(2)*vn*H0+aL(3)*vt*H0+0.5*aL(4)*((vn**2+vt**2)*H0+B0))*delta-&
                tau*Mt(2)*fn*vn*(aR(1)*H0+aR(2)*vn*H0+aR(3)*vt*H0+0.5*aR(4)*((vn**2+vt**2)*H0+B0))*(1-delta)+&
                Mt(6)*fn*vn*shn
    face%flux_b=face%flux_b-&
                tau*Mt(2)*fn*vn*(aL(1)*B0+aL(2)*vn*B0+aL(3)*vt*B0+0.5*aL(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*delta-&
                tau*Mt(2)*fn*vn*(aR(1)*B0+aR(2)*vn*B0+aR(3)*vt*B0+0.5*aR(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*(1-delta)+&
                Mt(6)*fn*vn*sbn

    face%flux_h=face%flux_h-&
                tau*Mt(2)*ft*vn*(ab(1)*H0+ab(2)*vn*H0+ab(3)*vt*H0+0.5*ab(4)*((vn**2+vt**2)*H0+B0))+&
                Mt(6)*ft*vn*sht
    face%flux_b=face%flux_b-&
                tau*Mt(2)*ft*vn*(ab(1)*B0+ab(2)*vn*B0+ab(3)*vt*B0+0.5*ab(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                Mt(6)*ft*vn*sbt

    !--------------------------------------------------
    ! final flux
    !--------------------------------------------------
    ! convert to global frame
    face%flux=global_frame(face%flux,face%cosa,face%sina) 
    ! total flux (乘以长度)
    face%flux=face%length*face%flux 
    face%flux_h=face%length*face%flux_h
    face%flux_b=face%length*face%flux_b

    end subroutine calc_flux_mtforcev
! ------------------------------------------------------------
    subroutine calc_flux_mtforceh(cell_L,face,cell_R,dirc)
    type(cell_interface),intent(inout):: face
    integer,intent(in):: dirc
    type(cell_center),intent(in):: cell_L,cell_R

    !Heaviside step function
    integer:: delta(unum,vnum)

    !Local micro velocity space
    real(kind=8):: vn(unum,vnum),vt(unum,vnum)

    !interface variable
    real(kind=8):: h(unum,vnum),b(unum,vnum)
    real(kind=8):: H0(unum,vnum),B0(unum,vnum)
    real(kind=8):: H_plus(unum,vnum),B_plus(unum,vnum)
    real(kind=8):: shn(unum,vnum),sbn(unum,vnum),sht(unum,vnum),sbt(unum,vnum)
    real(kind=8):: w(4),prim(4)
    real(kind=8):: qf(2)
    real(kind=8):: sw(4)
    real(kind=8):: aL(4),aR(4),aT(4),ab(4)

    !moments variable
    real(kind=8):: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
    real(kind=8):: Mau_0(4) !<u\psi>
    real(kind=8):: Mau_L(4),Mau_R(4) !<aL*u^n*\psi>,<aR*u^n*\psi>
    real(kind=8):: Mbu(4) !<b*u^n*\psi>
    real(kind=8):: Mau_T(4) !<A*u*\psi>
    real(kind=8):: tau
    real(kind=8):: Mt(6)

    !trace local velocity space
    real(kind=8):: vf(unum,vnum)
    real(kind=8):: fn,ft

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !--------------------------------------------------
    ! local variable
    !--------------------------------------------------
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina

    call face_acceleration(face,fn,ft)
    vf=vn-0.5d0*fn*dt

    delta=(sign(UP,vf)+1.d0)/2.d0

    !--------------------------------------------------
    ! reconstruct initial distribution
    !--------------------------------------------------
    ! UPWIND reconstruction
    h=(cell_L%h+0.5*cell_L%length(DIRC)*cell_L%sh(:,:,DIRC))*delta+&
      (cell_R%h-0.5*cell_R%length(DIRC)*cell_R%sh(:,:,DIRC))*(1-delta)
    b=(cell_L%b+0.5*cell_L%length(DIRC)*cell_L%sb(:,:,DIRC))*delta+&
      (cell_R%b-0.5*cell_R%length(DIRC)*cell_R%sb(:,:,DIRC))*(1-delta)
    shn=cell_L%sh(:,:,DIRC)*delta+cell_R%sh(:,:,DIRC)*(1-delta)
    ! For a horizontal face vt=-u, so the local tangential coordinate points
    ! in the negative global-x direction.
    sht=-cell_L%sh(:,:,mod(DIRC,2)+1)*delta-cell_R%sh(:,:,mod(DIRC,2)+1)*(1-delta)
    sbn=cell_L%sb(:,:,DIRC)*delta+cell_R%sb(:,:,DIRC)*(1-delta)
    sbt=-cell_L%sb(:,:,mod(DIRC,2)+1)*delta-cell_R%sb(:,:,mod(DIRC,2)+1)*(1-delta)

    !--------------------------------------------------
    ! obtain macroscopic variables at interface
    !--------------------------------------------------
    !conservative variables W_0 
    w(1)=sum(weight*h)
    w(2)=sum(weight*vn*h)
    w(3)=sum(weight*vt*h)
    w(4)=0.5*(sum(weight*(vn**2+vt**2)*h)+sum(weight*b))

    !convert to primary variables
    prim=get_primary(w)

    !heat flux
    qf=get_heat_flux(h,b,vn,vt,prim) 

    !--------------------------------------------------
    ! calculate a^L,a^R
    !--------------------------------------------------
    sw=(w-local_frame(cell_L%w,face%cosa,face%sina))/(0.5*cell_L%length(DIRC))
    aL=micro_slope(prim,sw)

    sw=(local_frame(cell_R%w,face%cosa,face%sina)-w)/(0.5*cell_R%length(DIRC))
    aR=micro_slope(prim,sw)

    !--------------------------------------------------
    ! calculate b
    !--------------------------------------------------
    sw(1)=sum(weight*sht)
    sw(2)=sum(weight*vn*sht)
    sw(3)=sum(weight*vt*sht)
    sw(4)=0.5d0*(sum(weight*(vn**2+vt**2)*sht)+sum(weight*sbt))
    ab=micro_slope(prim,sw)

    !--------------------------------------------------
    ! calculate time slope of W and A
    !--------------------------------------------------
    !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
    call calc_moment(prim,Mu,Mv,Mxi,Mu_L,Mu_R) 

    Mau_L=moment_au(aL,Mu_L,Mv,Mxi,1,0) !<aL*u*\psi>_{>0}
    Mau_R=moment_au(aR,Mu_R,Mv,Mxi,1,0) !<aR*u*\psi>_{<0}
    Mbu = moment_au(ab,Mu,Mv,Mxi,0,1)

    sw=-prim(1)*(Mau_L+Mau_R+Mbu) !time slope of W
    call reduced_maxwell(H0,B0,vn,vt,prim)
    call add_force_compatibility(sw,prim,vn,vt,H0,B0,fn,ft)
    aT=micro_slope(prim,sw) !calculate A

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
    Mau_0=moment_uv(Mu,Mv,Mxi,1,0,0) !<u*\psi>
    Mau_L=moment_au(aL,Mu_L,Mv,Mxi,2,0) !<aL*u^2*\psi>_{>0}
    Mau_R=moment_au(aR,Mu_R,Mv,Mxi,2,0) !<aR*u^2*\psi>_{<0}
    Mbu=moment_au(ab,Mu,Mv,Mxi,1,1) !!<b*u*v*\psi>
    Mau_T=moment_au(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>

    face%flux=Mt(1)*prim(1)*Mau_0+Mt(2)*prim(1)*(Mau_L+Mau_R+Mbu)+Mt(3)*prim(1)*Mau_T

    !--------------------------------------------------
    ! gravity flux related to g0
    !--------------------------------------------------
    !normal direction
    Mau_L=moment_au(aL,Mu_L,Mv,Mxi,1,0)
    Mau_R=moment_au(aR,Mu_R,Mv,Mxi,1,0)

    face%flux=face%flux-tau*Mt(2)*prim(1)*fn*(Mau_L+Mau_R)-&
              Mt(6)*prim(1)*fn*(Mau_L+Mau_R)
    !tangential direction
    Mbu=moment_au(ab,Mu,Mv,Mxi,0,1)
    face%flux=face%flux-tau*Mt(2)*prim(1)*ft*Mbu-&
              Mt(6)*prim(1)*ft*Mbu

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    !Maxwellian distribution H0 and B0
    call reduced_maxwell(H0,B0,vn,vt,prim)

    !Shakhov part H+ and B+
    call shakhov_part(H0,B0,vn,vt,qf,prim,H_plus,B_plus)

    !macro flux related to g+ and f0
    face%flux(1)=face%flux(1)+Mt(1)*sum(weight*vn*H_plus)+Mt(4)*sum(weight*vn*h)-Mt(5)*sum(weight*vn**2*shn)-Mt(5)*sum(weight*vt*vn*sht)
    face%flux(2)=face%flux(2)+Mt(1)*sum(weight*vn*vn*H_plus)+Mt(4)*sum(weight*vn*vn*h)-Mt(5)*sum(weight*vn*vn**2*shn)-Mt(5)*sum(weight*vn*vt*vn*sht)
    face%flux(3)=face%flux(3)+Mt(1)*sum(weight*vt*vn*H_plus)+Mt(4)*sum(weight*vt*vn*h)-Mt(5)*sum(weight*vt*vn**2*shn)-Mt(5)*sum(weight*vt*vt*vn*sht)
    face%flux(4)=face%flux(4)+&
                 Mt(1)*0.5*(sum(weight*vn*(vn**2+vt**2)*H_plus)+sum(weight*vn*B_plus))+&
                 Mt(4)*0.5*(sum(weight*vn*(vn**2+vt**2)*h)+sum(weight*vn*b))-&
                 Mt(5)*0.5*(sum(weight*vn**2*(vn**2+vt**2)*shn)+sum(weight*vn**2*sbn))-&
                 Mt(5)*0.5*(sum(weight*vt*vn*(vn**2+vt**2)*sht)+sum(weight*vt*vn*sbt))

    !gravity related to g+ and f0
    face%flux(1)=face%flux(1)+fn*Mt(6)*sum(weight*vn*shn)
    face%flux(2)=face%flux(2)+fn*Mt(6)*sum(weight*vn*vn*shn)
    face%flux(3)=face%flux(3)+fn*Mt(6)*sum(weight*vt*vn*shn)
    face%flux(4)=face%flux(4)+fn*Mt(6)*0.5*(sum(weight*(vn**2+vt**2)*vn*shn)+sum(weight*vn*sbn))

    face%flux(1)=face%flux(1)+ft*Mt(6)*sum(weight*vn*sht)
    face%flux(2)=face%flux(2)+ft*Mt(6)*sum(weight*vn*vn*sht)
    face%flux(3)=face%flux(3)+ft*Mt(6)*sum(weight*vt*vn*sht)
    face%flux(4)=face%flux(4)+ft*Mt(6)*0.5*(sum(weight*(vn**2+vt**2)*vn*sht)+sum(weight*vn*sbt))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    face%flux_h=Mt(1)*vn*(H0+H_plus)+&
                Mt(2)*vn**2*(aL(1)*H0+aL(2)*vn*H0+aL(3)*vt*H0+0.5*aL(4)*((vn**2+vt**2)*H0+B0))*delta+& 
                Mt(2)*vn**2*(aR(1)*H0+aR(2)*vn*H0+aR(3)*vt*H0+0.5*aR(4)*((vn**2+vt**2)*H0+B0))*(1-delta)+&
                Mt(2)*vn*vt*(ab(1)*H0+ab(2)*vn*H0+ab(3)*vt*H0+0.5*ab(4)*((vn**2+vt**2)*H0+B0))+&
                Mt(3)*vn*(aT(1)*H0+aT(2)*vn*H0+aT(3)*vt*H0+0.5*aT(4)*((vn**2+vt**2)*H0+B0))+&
                Mt(4)*vn*h-Mt(5)*vn**2*shn-Mt(5)*vn*vt*sht

    face%flux_b=Mt(1)*vn*(B0+B_plus)+&
                Mt(2)*vn**2*(aL(1)*B0+aL(2)*vn*B0+aL(3)*vt*B0+0.5*aL(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*delta+&
                Mt(2)*vn**2*(aR(1)*B0+aR(2)*vn*B0+aR(3)*vt*B0+0.5*aR(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*(1-delta)+&
                Mt(2)*vn*vt*(ab(1)*B0+ab(2)*vn*B0+ab(3)*vt*B0+0.5*ab(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                Mt(3)*vn*(aT(1)*B0+aT(2)*vn*B0+aT(3)*vt*B0+0.5*aT(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                Mt(4)*vn*b-Mt(5)*vn**2*sbn-Mt(5)*vn*vt*sbt

    !gravity for distribution function
    face%flux_h=face%flux_h-&
                tau*Mt(2)*fn*vn*(aL(1)*H0+aL(2)*vn*H0+aL(3)*vt*H0+0.5*aL(4)*((vn**2+vt**2)*H0+B0))*delta-&
                tau*Mt(2)*fn*vn*(aR(1)*H0+aR(2)*vn*H0+aR(3)*vt*H0+0.5*aR(4)*((vn**2+vt**2)*H0+B0))*(1-delta)+&
                Mt(6)*fn*vn*shn
    face%flux_b=face%flux_b-&
                tau*Mt(2)*fn*vn*(aL(1)*B0+aL(2)*vn*B0+aL(3)*vt*B0+0.5*aL(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*delta-&
                tau*Mt(2)*fn*vn*(aR(1)*B0+aR(2)*vn*B0+aR(3)*vt*B0+0.5*aR(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*(1-delta)+&
                Mt(6)*fn*vn*sbn

    face%flux_h=face%flux_h-&
                tau*Mt(2)*ft*vn*(ab(1)*H0+ab(2)*vn*H0+ab(3)*vt*H0+0.5*ab(4)*((vn**2+vt**2)*H0+B0))+&
                Mt(6)*ft*vn*sht
    face%flux_b=face%flux_b-&
                tau*Mt(2)*ft*vn*(ab(1)*B0+ab(2)*vn*B0+ab(3)*vt*B0+0.5*ab(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                Mt(6)*ft*vn*sbt

    !--------------------------------------------------
    ! final flux
    !--------------------------------------------------
    ! convert to global frame
    face%flux=global_frame(face%flux,face%cosa,face%sina) 
    ! total flux (乘以长度)
    face%flux=face%length*face%flux 
    face%flux_h=face%length*face%flux_h
    face%flux_b=face%length*face%flux_b

    end subroutine calc_flux_mtforceh
! ------------------------------------------------------------
    subroutine calc_flux_boundary(bc,face,cell,DIRC,rot) 
    real(kind=8),intent(in):: bc(DIM+2)
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell
    integer,intent(in):: DIRC,rot
    
    !Heaviside step function
    integer:: delta(unum,vnum)
    
    !local micro velocity space
    real(kind=8):: vn(unum,vnum),vt(unum,vnum)
    
    !interface variable
    real(kind=8):: h(unum,vnum),b(unum,vnum)
    real(kind=8):: H0(unum,vnum),B0(unum,vnum)
    real(kind=8):: prim(4) !boundary condition in local frame
    real(kind=8):: R1,R2 !±ﬂΩÁ√‹∂»º∆À„¡ø

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !allocate array
    !convert the micro velocity to local frame
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina

    !Heaviside step function. The rotation accounts for the right wall
    delta=(sign(UP,vn)*rot+1)/2

    !boundary condition in local frame
    prim=local_frame(bc,face%cosa,face%sina)

    !--------------------------------------------------
    !obtain h^{in} and b^{in}, rotation accounts for the right wall
    !--------------------------------------------------
    h=cell%h-rot*0.5*cell%length(DIRC)*cell%sh(:,:,DIRC)
    b=cell%b-rot*0.5*cell%length(DIRC)*cell%sb(:,:,DIRC)

    !--------------------------------------------------
    !calculate wall density and Maxwellian distribution
    !--------------------------------------------------
    R1=sum(weight*vn*h*(1-delta))
    R2=(prim(4)/PI)*sum(weight*vn*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))*delta)

    prim(1)=-R1/R2

    call reduced_maxwell(H0,B0,vn,vt,prim)

    !--------------------------------------------------
    !distribution function at the boundary interface
    !--------------------------------------------------
    h=H0*delta+h*(1-delta)
    b=B0*delta+b*(1-delta)

    !--------------------------------------------------
    !calculate flux
    !--------------------------------------------------
    face%flux(1)=sum(weight*vn*h)
    face%flux(2)=sum(weight*vn*vn*h)
    face%flux(3)=sum(weight*vn*vt*h)
    face%flux(4)=0.5*sum(weight*vn*((vn**2+vt**2)*h+b))

    face%flux_h=vn*h
    face%flux_b=vn*b

    !--------------------------------------------------
    !final flux
    !--------------------------------------------------
    !convert to global frame
    face%flux=global_frame(face%flux,face%cosa,face%sina) 
    !total flux
    face%flux=dt*face%length*face%flux
    face%flux_h=dt*face%length*face%flux_h
    face%flux_b=dt*face%length*face%flux_b

    end subroutine calc_flux_boundary
! ------------------------------------------------------------
    subroutine calc_flux_specular_vertical(face,cell,DIRC,rot) 
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell
    integer,intent(in):: DIRC,rot
    
    !Heaviside step function
    integer:: delta(unum,vnum)
    
    !local micro velocity space
    real(kind=8):: vn(unum,vnum),vt(unum,vnum)
    
    !interface variable
    real(kind=8):: h(unum,vnum),b(unum,vnum)
    real(kind=8):: prim(4) !boundary condition in local frame
    real(kind=8):: temph(unum,vnum),tempb(unum,vnum)

    integer:: i,j

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !allocate array
    !convert the micro velocity to local frame
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina

    !Heaviside step function. The rotation accounts for the right wall
    delta=(sign(UP,vn)*rot+1)/2

    !--------------------------------------------------
    !obtain h^{in} and b^{in}, rotation accounts for the right wall
    !--------------------------------------------------
    h=cell%h-rot*0.5*cell%length(DIRC)*cell%sh(:,:,DIRC)
    b=cell%b-rot*0.5*cell%length(DIRC)*cell%sb(:,:,DIRC)

    !--------------------------------------------------
    !specular distribution function at the boundary interface
    !--------------------------------------------------
    do i=1,unum
    do j=1,vnum
        temph(i,j)=h(unum-i+1,j)
        tempb(i,j)=b(unum-i+1,j)
    end do
    end do
    h=temph*delta+h*(1-delta)
    b=tempb*delta+b*(1-delta)

    !--------------------------------------------------
    !calculate flux
    !--------------------------------------------------
    face%flux(1)=sum(weight*vn*h)
    face%flux(2)=sum(weight*vn*vn*h)
    face%flux(3)=sum(weight*vn*vt*h)
    face%flux(4)=0.5*sum(weight*vn*((vn**2+vt**2)*h+b))

    face%flux_h=vn*h
    face%flux_b=vn*b

    !--------------------------------------------------
    !final flux
    !--------------------------------------------------
    !convert to global frame
    face%flux=global_frame(face%flux,face%cosa,face%sina) 
    !total flux
    face%flux=dt*face%length*face%flux
    face%flux_h=dt*face%length*face%flux_h
    face%flux_b=dt*face%length*face%flux_b

    end subroutine calc_flux_specular_vertical
! ------------------------------------------------------------
    subroutine calc_flux_specular_horizon(face,cell,DIRC,rot) 
    type(cell_interface),intent(inout):: face
    type(cell_center),intent(in):: cell
    integer,intent(in):: DIRC,rot
    
    !Heaviside step function
    integer:: delta(unum,vnum)
    
    !local micro velocity space
    real(kind=8):: vn(unum,vnum),vt(unum,vnum)
    
    !interface variable
    real(kind=8):: h(unum,vnum),b(unum,vnum)
    real(kind=8):: prim(4) !boundary condition in local frame
    real(kind=8):: temph(unum,vnum),tempb(unum,vnum)

    integer:: i,j

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    !allocate array
    !convert the micro velocity to local frame
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina

    !Heaviside step function. The rotation accounts for the right wall
    delta=(sign(UP,vn)*rot+1)/2

    !--------------------------------------------------
    !obtain h^{in} and b^{in}, rotation accounts for the right wall
    !--------------------------------------------------
    h=cell%h-rot*0.5*cell%length(DIRC)*cell%sh(:,:,DIRC)
    b=cell%b-rot*0.5*cell%length(DIRC)*cell%sb(:,:,DIRC)

    !--------------------------------------------------
    !specular distribution function at the boundary interface
    !--------------------------------------------------
    do i=1,unum
    do j=1,vnum
        temph(i,j)=h(i,vnum-j+1)
        tempb(i,j)=b(i,vnum-j+1)
    end do
    end do
    h=temph*delta+h*(1-delta)
    b=tempb*delta+b*(1-delta)

    !--------------------------------------------------
    !calculate flux
    !--------------------------------------------------
    face%flux(1)=sum(weight*vn*h)
    face%flux(2)=sum(weight*vn*vn*h)
    face%flux(3)=sum(weight*vn*vt*h)
    face%flux(4)=0.5*sum(weight*vn*((vn**2+vt**2)*h+b))

    face%flux_h=vn*h
    face%flux_b=vn*b

    !--------------------------------------------------
    !final flux
    !--------------------------------------------------
    !convert to global frame
    face%flux=global_frame(face%flux,face%cosa,face%sina) 
    !total flux
    face%flux=dt*face%length*face%flux
    face%flux_h=dt*face%length*face%flux_h
    face%flux_b=dt*face%length*face%flux_b

    end subroutine calc_flux_specular_horizon
! ------------------------------------------------------------
    subroutine add_force_compatibility(sw,prim,vn,vt,H0,B0,fn,ft)
        ! Compatibility condition for the time derivative of the
        ! equilibrium distribution in an external acceleration field.
        real(kind=8),intent(inout):: sw(DIM+2)
        real(kind=8),intent(in):: prim(DIM+2),vn(:,:),vt(:,:),H0(:,:),B0(:,:)
        real(kind=8),intent(in):: fn,ft

        sw(2)=sw(2)+2.d0*sum(weight*vn*prim(4)* &
              (fn*(vn-prim(2))+ft*(vt-prim(3)))*H0)
        sw(3)=sw(3)+2.d0*sum(weight*vt*prim(4)* &
              (fn*(vn-prim(2))+ft*(vt-prim(3)))*H0)
        sw(4)=sw(4)+sum(weight*prim(4)* &
              (fn*(vn-prim(2))+ft*(vt-prim(3)))* &
              ((vn**2+vt**2)*H0+B0))
    end subroutine add_force_compatibility
! ------------------------------------------------------------
    !--------------------------------------------------
    ! calculate slope of Maxwellian distribution
    ! sw(4): slope of macro variables W
    !--------------------------------------------------
    function micro_slope(prim,sw)
        real(kind=8),intent(in):: prim(DIM+2),sw(DIM+2)
        real(kind=8):: micro_slope(DIM+2)

        micro_slope(4)=4.0*prim(4)**2/(innerk+2)/prim(1)*(2.0*sw(4)-2.0*prim(2)*sw(2)-&
        		   2.0*prim(3)*sw(3)+sw(1)*(prim(2)**2+prim(3)**2-0.5*(innerk+2)/prim(4)))

        micro_slope(3)=2.0*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*micro_slope(4)
        micro_slope(2)=2.0*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*micro_slope(4)
        micro_slope(1)=sw(1)/prim(1)-prim(2)*micro_slope(2)-prim(3)*micro_slope(3)-&
        		   0.5*(prim(2)**2+prim(3)**2+0.5*(innerk+2)/prim(4))*micro_slope(4)
    end function micro_slope
! ------------------------------------------------------------
    !--------------------------------------------------
    ! calculate moments of micro velocity and internal freedom 
    ! Mu,Mv     :<u^n>,<v^m>
    ! Mxi       :<\xi^l>
    ! Mu_L,Mu_R :<u^n>_{>0},<u^n>_{<0}
    !--------------------------------------------------
    subroutine calc_moment(prim,Mu,Mv,Mxi,Mu_L,Mu_R)
        real(kind=8),intent(in):: prim(DIM+2)
        real(kind=8),intent(out):: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)
        real(kind=8),intent(out):: Mv(0:MTUM)
        real(kind=8),intent(out):: Mxi(0:2)
        integer:: i

        !moments of normal velocity
        Mu_L(0)=0.5*erfc(-sqrt(prim(4))*prim(2))
        Mu_L(1)=prim(2)*Mu_L(0)+0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))
        Mu_R(0)=0.5*erfc(sqrt(prim(4))*prim(2))
        Mu_R(1)=prim(2)*Mu_R(0)-0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))

        do i=2,MNUM
            Mu_L(i)=prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)/prim(4)
            Mu_R(i)=prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)/prim(4)
        end do

        Mu=Mu_L+Mu_R

        !moments of tangential velocity
        Mv(0)=1.0
        Mv(1)=prim(3)

        do i=2,MTUM
            Mv(i)=prim(3)*Mv(i-1)+0.5*(i-1)*Mv(i-2)/prim(4)
        end do

        !moments of \xi
        Mxi(0)=1.0 !<\xi^0>
        Mxi(1)=0.5*innerk/prim(4) !<\xi^2>
        Mxi(2)=(innerk**2+2.0*innerk)/(4.0*prim(4)**2) !<\xi^4>
    end subroutine calc_moment
! ------------------------------------------------------------
    !--------------------------------------------------
    ! calculate <u^\alpha*v^\beta*\xi^\theta*\psi>
    ! Mu,Mv      :<u^\alpha>,<v^\beta>
    ! Mxi        :<\xi^l>
    ! alpha,beta :exponential index of u and v
    ! theta      :exponential index of \xi
    !moment_uv   :moment of <u^\alpha*v^\beta*\xi^\theta*\psi>
    !--------------------------------------------------
    function moment_uv(Mu,Mv,Mxi,alpha,beta,theta)
        real(kind=8),intent(in):: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer,intent(in):: alpha,beta,theta
        real(kind=8):: moment_uv(DIM+2)

        moment_uv(1)=Mu(alpha)*Mv(beta)*Mxi(theta/2)
        moment_uv(2)=Mu(alpha+1)*Mv(beta)*Mxi(theta/2)
        moment_uv(3)=Mu(alpha)*Mv(beta+1)*Mxi(theta/2)
        moment_uv(4)=0.5*(Mu(alpha+2)*Mv(beta)*Mxi(theta/2)+&
        		 Mu(alpha)*Mv(beta+2)*Mxi(theta/2)+Mu(alpha)*Mv(beta)*Mxi((theta+2)/2))
    end function moment_uv
! ------------------------------------------------------------
    !--------------------------------------------------
    ! calculate <a*u^\alpha*v^\beta*\psi>
    ! a          :micro slope of Maxwellian
    ! Mu,Mv      :<u^\alpha>,<v^\beta>
    ! Mxi        :<\xi^l>
    ! alpha,beta :exponential index of u and v
    ! moment_au  :moment of <a*u^\alpha*v^\beta*\psi>
    !--------------------------------------------------
    function moment_au(a,Mu,Mv,Mxi,alpha,beta)
        real(kind=8),intent(in):: a(DIM+2)
        real(kind=8),intent(in):: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer,intent(in):: alpha,beta
        real(kind=8):: moment_au(DIM+2)

        moment_au = a(1)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,0)+&
                    a(2)*moment_uv(Mu,Mv,Mxi,alpha+1,beta+0,0)+&
                    a(3)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+1,0)+&
                    0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+2,beta+0,0)+&
                    0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+2,0)+&
                    0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,2)
    end function moment_au
! ------------------------------------------------------------
end module FLUX
