! ------------------------------------------------------------
! Gravity source term in the distribution function update
! CONTINUOUS distribution function --> UPWIND/CENTRAL difference
! vanLeer limiter 
! ------------------------------------------------------------
module GSOURCE
use ALLDATA
use GASTHEORY
implicit none

real(kind=8),pointer:: duh(:,:)
real(kind=8),pointer:: dub(:,:)
real(kind=8),pointer:: dvh(:,:)
real(kind=8),pointer:: dvb(:,:)


contains
! ------------------------------------------------------------
	subroutine cal_fsource(cell,face)
    type(cell_interface),intent(inout):: face
	type(cell_center),intent(in):: cell
	integer:: i,j
    real(kind=8),allocatable,dimension(:,:):: vn,vt
	real(kind=8):: sL,sR	
	real(kind=8):: du,dv
    real(kind=8):: w1(4),prim1(4)
	if(associated(duh)) deallocate(duh)
	allocate(duh(unum,vnum))
	if(associated(dub)) deallocate(dub)
	allocate(dub(unum,vnum))
	if(associated(dvh)) deallocate(dvh)
	allocate(dvh(unum,vnum))
	if(associated(dvb)) deallocate(dvb)
	allocate(dvb(unum,vnum))
    allocate(vn(unum,vnum))
    allocate(vt(unum,vnum))
    
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina
    
    w1(1)=sum(weight*(cell%h))
    w1(2)=sum(weight*vn*(cell%h))
    w1(3)=sum(weight*vt*(cell%h))
    w1(4)=0.5*(sum(weight*(vn**2+vt**2)*(cell%h))+sum(weight*(cell%b)))
    prim1=get_primary(w1)

	du=lengthU/unum
	dv=lengthV/vnum

	!boundary
	duh(1,:)=(cell%h(2,:)-cell%h(1,:))/du
	dub(1,:)=(cell%b(2,:)-cell%b(1,:))/du
	dvh(:,1)=(cell%h(:,2)-cell%h(:,1))/dv
	dvb(:,1)=(cell%b(:,2)-cell%b(:,1))/dv
	duh(unum,:)=(cell%h(unum,:)-cell%h(unum-1,:))/du
	dub(unum,:)=(cell%b(unum,:)-cell%b(unum-1,:))/du
	dvh(:,vnum)=(cell%h(:,vnum)-cell%h(:,vnum-1))/dv
	dvb(:,vnum)=(cell%b(:,vnum)-cell%b(:,vnum-1))/dv

	!inner
	do i=2,unum-1
	do j=1,vnum
			sl=(cell%h(i,j)-cell%h(i-1,j))/du
			sr=(cell%h(i+1,j)-cell%h(i,j))/du
			duh(i,j)=sr
			!duh(i,j)=0.5*(sl+sr)
			!duh(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)

			sl=(cell%b(i,j)-cell%b(i-1,j))/du
			sr=(cell%b(i+1,j)-cell%b(i,j))/du
			dub(i,j)=sr
			!dub(i,j)=0.5*(sl+sr)
			!dub(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)
          !duh(i,j)=-2*(vn(i,j)-prim1(2))*prim1(4)*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
          !dub(i,j)=-1*(vn(i,j)-prim1(2))*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
	end do
	end do

	do i=1,unum
	do j=2,vnum-1
			sl=(cell%h(i,j)-cell%h(i,j-1))/dv
			sr=(cell%h(i,j+1)-cell%h(i,j))/dv
			dvh(i,j)=sr
			!dvh(i,j)=0.5*(sl+sr)
			!dvh(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)
			
			sl=(cell%b(i,j)-cell%b(i,j-1))/dv
			sr=(cell%b(i,j+1)-cell%b(i,j))/dv
			dvb(i,j)=sr
			!dvb(i,j)=0.5*(sl+sr)
			!dvb(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)
            !dvh(i,j)=-2*(vt(i,j)-prim1(3))*prim1(4)*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
            !dvb(i,j)=-1*(vt(i,j)-prim1(3))*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
	end do
	end do

    end subroutine cal_fsource
! ------------------------------------------------------
subroutine cal_fsource1(cell,face)
    type(cell_interface),intent(inout):: face
	type(cell_center),intent(in):: cell
	integer:: i,j
    real(kind=8),allocatable,dimension(:,:):: vn,vt
	real(kind=8):: sL,sR	
	real(kind=8):: du,dv
    real(kind=8):: w1(4),prim1(4)
	if(associated(duh)) deallocate(duh)
	allocate(duh(unum,vnum))
	if(associated(dub)) deallocate(dub)
	allocate(dub(unum,vnum))
	if(associated(dvh)) deallocate(dvh)
	allocate(dvh(unum,vnum))
	if(associated(dvb)) deallocate(dvb)
	allocate(dvb(unum,vnum))
    allocate(vn(unum,vnum))
    allocate(vt(unum,vnum))
    
    vn=uspace*face%cosa+vspace*face%sina
    vt=vspace*face%cosa-uspace*face%sina
    
    w1(1)=sum(weight*(cell%h))
    w1(2)=sum(weight*vn*(cell%h))
    w1(3)=sum(weight*vt*(cell%h))
    w1(4)=0.5*(sum(weight*(vn**2+vt**2)*(cell%h))+sum(weight*(cell%b)))
    prim1=get_primary(w1)

	du=lengthU/unum
	dv=lengthV/vnum

	!boundary
	duh(1,:)=(-3.0d0*cell%h(1,:)+4.0d0*cell%h(2,:)-cell%h(3,:))/du/2.0d0
    dub(1,:)=(-3.0d0*cell%b(1,:)+4.0d0*cell%b(2,:)-cell%b(3,:))/du/2.0d0
    duh(2,:)=(cell%h(3,:)-cell%h(1,:))/du/2.0d0
    dub(2,:)=(cell%b(3,:)-cell%b(1,:))/du/2.0d0
    dvh(:,1)=(-3.0d0*cell%h(:,1)+4.0d0*cell%h(:,2)-cell%h(:,3))/dv/2.0d0
    dvb(:,1)=(-3.0d0*cell%b(:,1)+4.0d0*cell%b(:,2)-cell%b(:,3))/dv/2.0d0
    dvh(:,2)=(cell%h(:,3)-cell%h(:,1))/dv/2.0d0
    dvb(:,2)=(cell%b(:,3)-cell%b(:,1))/dv/2.0d0
    
    
	duh(unum,:)=(3.0d0*cell%h(unum,:)-4.0d0*cell%h(unum-1,:)+cell%h(unum-2,:))/du/2.0d0
    dub(unum,:)=(3.0d0*cell%b(unum,:)-4.0d0*cell%b(unum-1,:)+cell%b(unum-2,:))/du/2.0d0
    duh(unum-1,:)=(cell%h(unum,:)-cell%h(unum-2,:))/du/2.0d0
    dub(unum-1,:)=(cell%b(unum,:)-cell%b(unum-2,:))/du/2.0d0
    dvh(:,unum)=(-3.0d0*cell%h(:,unum)+4.0d0*cell%h(:,unum-1)-cell%h(:,unum-2))/dv/2.0d0
    dvb(:,unum)=(-3.0d0*cell%b(:,unum)+4.0d0*cell%b(:,unum-1)-cell%b(:,unum-2))/dv/2.0d0
    dvh(:,unum-1)=(cell%h(:,unum)-cell%h(:,unum-2))/dv/2.0d0
    dvb(:,unum-1)=(cell%b(:,unum)-cell%b(:,unum-2))/dv/2.0d0

	!inner
	do i=3,unum-2
	do j=1,vnum
			sl=(cell%h(i,j)-cell%h(i-1,j))/du
			sr=(-cell%h(i+2,j)+8.0d0*cell%h(i+1,j)-8.0d0*cell%h(i-1,j)+cell%h(i-2,j))/du/12.0d0
			duh(i,j)=sr
			!duh(i,j)=0.5*(sl+sr)
			!duh(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)

			sl=(cell%b(i,j)-cell%b(i-1,j))/du
			sr=(-cell%b(i+2,j)+8.0d0*cell%h(i+1,j)-8.0d0*cell%h(i-1,j)+cell%h(i-2,j))/du/12.0d0
			dub(i,j)=sr
			!dub(i,j)=0.5*(sl+sr)
			!dub(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)
          !duh(i,j)=-2*(vn(i,j)-prim1(2))*prim1(4)*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
          !dub(i,j)=-1*(vn(i,j)-prim1(2))*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
	end do
	end do

	do i=1,unum
	do j=3,vnum-2
			sl=(cell%h(i,j)-cell%h(i,j-1))/dv
			sr=(-cell%h(i,j+2)+8.0d0*cell%h(i,j+1)-8.0d0*cell%h(i,j-1)+cell%h(i,j-2))/dv/12.0d0
			dvh(i,j)=sr
			!dvh(i,j)=0.5*(sl+sr)
			!dvh(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)
			
			sl=(cell%b(i,j)-cell%b(i,j-1))/dv
			sr=(-cell%b(i,j+2)+8.0d0*cell%b(i,j+1)-8.0d0*cell%b(i,j-1)+cell%b(i,j-2))/dv/12.0d0
			dvb(i,j)=sr
			!dvb(i,j)=0.5*(sl+sr)
			!dvb(i,j)=(sign(up,sl)+sign(up,sr))*abs(sl)*abs(sr)/(abs(sl)+abs(sr)+smv)
            !dvh(i,j)=-2*(vt(i,j)-prim1(3))*prim1(4)*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
            !dvb(i,j)=-1*(vt(i,j)-prim1(3))*prim1(1)*(prim1(4)/PI)*EXP(-prim1(4)*((vn(i,j)-prim1(2))**2+(vt(i,j)-prim1(3))**2))
	end do
	end do

end subroutine cal_fsource1
! ------------------------------------------------------
subroutine new_expansion_force_term(h,b,order,prim , g_x , g_y, force_term, force_term_t,error,errort)
    implicit none
    integer, INTENT(in) :: order
    real(kind=8):: prim(DIM+2)
    REAL*8, INTENT(in):: h(unum, vnum), b(unum,vnum), g_x, g_y
    real*8, intent(out) :: error, errort,force_term(unum,vnum), force_term_t(unum,vnum)
    REAL*8 :: w(unum,vnum), h_ex(unum,vnum)
    REAL*8 :: vx2d(unum,vnum),vy2d(unum,vnum)
    REAL*8 :: dn(0:order,0:order)

    real*8 :: frac(0:order), vx(unum), vy(vnum), theta
    integer :: i, j
    real*8 :: hpx(0:order+1,unum,vnum),hpy(0:order+1,unum,vnum)

    theta=1.0d0/prim(4)
    
    vx2d = (uspace-prim(2))/dsqrt(theta)*dsqrt(2.d0)
    vy2d = (vspace-prim(3))/dsqrt(theta)*dsqrt(2.d0)

    w = 1.d0/2.d0/pi*dexp(-( vx2d**2+ vy2d**2)/2.d0)

    frac(0) = 1.d0
    do I = 1,order
        frac(i) = dble(i)*frac(i-1)     
    end do

    hpx(0,:,:) = 1.d0
    hpy(0,:,:) = 1.d0
    hpx(1,:,:) = vx2d
    hpy(1,:,:) = vy2d
    if ( order.ge.1 ) then
        do I = 2, order+1
            hpx(i,:,:) = vx2d*hpx(i-1,:,:)-(i-1.d0)*hpx(i-2,:,:)
            hpy(i,:,:) = vy2d*hpy(i-1,:,:)-(i-1.d0)*hpy(i-2,:,:)
        end do
    end if


    dn = 0.d0
    dn(0,0) = sum(h*hpx(0,:,:)*hpy(0,:,:)*weight) * 2.d0/theta  
    do I = 0, order
        do j = 0, i
            dn(i,j) = sum(h*hpx(j,:,:)*hpy(i-j,:,:)*weight) * 2.d0/theta  
        end do
    end do
    force_term = 0.d0
    h_ex = 0.d0
    do i = 0, order
        do j = 0, i
            force_term = force_term - dsqrt(2.d0/theta)*w*1.d0/frac(j)/frac(i-j)*&
            (g_x*hpx(j+1,:,:)*hpy(i-j,:,:)+g_y*hpx(j,:,:)*hpy(i-j+1,:,:))*dn(i,j)
            h_ex = h_ex + w*1.d0/frac(j)/frac(i-j)*hpx(j,:,:)*hpy(i-j,:,:)*dn(i,j)
        end do
    end do

    ! force_term = - dsqrt(2.d0/theta)*w*1.d0*&
    ! (g_x*hpx(1,:,:)*hpy(0,:,:)+g_y*hpx(0,:,:)*hpy(1,:,:))*dn(0,0)
    ! force_term = -1.d0/pi/theta*sum(h*weight)*dexp(-vx2d**2/2.d0-vy2d**2/2.d0)*(g_x*(uspace-u(1))+g_y*(vspace-u(2)))*2.d0/theta
    ! error = dsqrt(sum((h-h_ex)**2))/dsqrt(sum(h**2))
    error = dsqrt(sum((h-h_ex)**2))/dsqrt(sum(h**2))
    dn = 0.d0
    do I = 0, order
        do j = 0, i
            dn(i,j) = sum(b*hpx(j,:,:)*hpy(i-j,:,:)*weight) * 2.d0/theta  
        end do
    end do
    dn(0,0) = sum(b*hpx(0,:,:)*hpy(0,:,:)*weight) * 2.d0/theta  
    force_term_t = 0.d0
    h_ex = 0.d0
    do i = 0, order
        do j = 0, i
            force_term_t = force_term_t - dsqrt(2.d0/theta)*w*1.d0/frac(j)/frac(i-j)*&
            (g_x*hpx(j+1,:,:)*hpy(i-j,:,:)+g_y*hpx(j,:,:)*hpy(i-j+1,:,:))*dn(i,j)
            h_ex = h_ex + w*1.d0/frac(j)/frac(i-j)*hpx(j,:,:)*hpy(i-j,:,:)*dn(i,j)
        end do
    end do
    ! force_term_t = - dsqrt(2.d0/theta)*w*1.d0*&
    ! (g_x*hpx(1,:,:)*hpy(0,:,:)+g_y*hpx(0,:,:)*hpy(1,:,:))*dn(0,0)

    ! force_term_t = -1.d0/pi/2.d0*sum(b*weight)*dexp(-vx2d**2/2.d0-vy2d**2/2.d0)*(g_x*(uspace-u(1))+g_y*(vspace-u(2)))*2.d0/theta


    errort = dsqrt(sum((b-h_ex)**2))/dsqrt(sum(b**2))

end subroutine new_expansion_force_term
!-------------------------------------------------------------------
end module GSOURCE
