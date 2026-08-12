module GSOURCE
use CONFIG
use ALLDATA

implicit none

real(kind=8),pointer:: dvh(:)
real(kind=8),pointer:: dvb(:)
real(kind=8),pointer:: hH(:)
real(kind=8),pointer:: bH(:)
real(kind=8):: error1=0.d0
real(kind=8):: error2=0.d0
integer :: N_order,N_order1=0
public N_order1,error1,error2

contains
! ------------------------------------------------------------
	subroutine cal_fv(cell,prim)
	type(cell_center),intent(in):: cell
    real(kind=8):: prim(DIM+2)
    real(kind=8):: w1(DIM+2)
	real(kind=8):: du
    real(kind=8), allocatable :: factor(:)
    real(kind=8), allocatable :: hp(:,:)
    real(kind=8), allocatable :: F_term1(:,:)
    real(kind=8), allocatable :: F_term2(:,:)
    real(kind=8), allocatable :: v(:)
    real(kind=8), allocatable :: w(:)
    real(kind=8), allocatable :: df(:)
    real(kind=8), allocatable :: dg(:)
    
    real(kind=8):: rho,u,t
	integer:: i
    integer :: j
    
	if(associated(dvh)) deallocate(dvh)
	allocate(dvh(unum))
	if(associated(dvb)) deallocate(dvb)
	allocate(dvb(unum))
    if(associated(hH)) deallocate(hH)
	
    allocate( factor(0:N_order1))
    allocate( df(0:N_order1))
    allocate( dg(0:N_order1))
    allocate( hp(0:N_order1+1,unum))
    allocate( F_term1(unum,0:N_order1))
    allocate( F_term2(unum,0:N_order1))
    allocate(v(unum))
    allocate(w(unum))
 
    
    
    factor(0) = 1.0
    do i=1,N_order1
        factor(i) = 1.d0
        do j = 1 , i
            factor(i) = factor(i) * j
        end do
    end do
    
   
    rho = prim(1)
    u = prim(2)
    t = 1.0d0 / prim(3)
    v = (uspace-u)/dsqrt(t)*dsqrt(2.0d0)
    w = dexp(-v**2/2.d0)/dsqrt(2.d0*PI)
    !print *, "t = ",w1(1)
    error1=1.0d0
    error2=1.0d0
     
       hp(0,:) = 1.d0
       if ( N_order1 > -1 ) then
         hp(1,:) = v
       end if
      if ( N_order1 > 0 ) then
        hp(2,:) = v**2 - 1.d0
      end if
     if ( N_order1 > 1 ) then
       do j=3,N_order1+1
         hp(j,:)=v * hp(j-1,:) - ( j - 1.d0 ) * hp(j-2,:)
       end do
     end if
  ! do j=0,N_order
       ! do i=1,unum
       !     print *, "hp(", j, ",", i, ") = ", hp(1,1)
       ! end do
   !end do
  ! print *, "hp(", j, ",", i, ") = ", hp(1,1)
     do j=0,N_order1         
        df(j)=sum(cell%h * weight * hp(j,:))*dsqrt(2.d0/t)
        dg(j)=sum(cell%b * weight * hp(j,:))*dsqrt(2.d0/t)
         !print *, "df(", j, ") = ", df(j)
        !print *, "dg(", j, ") = ", dg(j)
     end do
    
     F_term1(:,:) = 0.d0
     F_term2(:,:) = 0.d0

     do i=1,unum
        do j=0,N_order1
            if (j < 1) then
                 F_term1(i,j) = F_term1(i,j)+df(j)*hp(j+1,i)/factor(j)
                 F_term2(i,j) = F_term2(i,j)+dg(j)*hp(j+1,i)/factor(j)
            end if
            if (j > 0) then
                 F_term1(i,j) = F_term1(i,j-1)+df(j)*hp(j+1,i)/factor(j)
                 F_term2(i,j) = F_term2(i,j-1)+dg(j)*hp(j+1,i)/factor(j)
        
            end if
           !  print *, "F_term1(", i, ",", j, ") = ", F_term1(i,j)
        end do
        dvh(i) =-dsqrt(2.0d0/t) * w(i) * F_term1(i,N_order1)
        dvb(i) =-dsqrt(2.0d0/t) * w(i) * F_term2(i,N_order1)
        ! dvh(i)=-2*(uspace(i)-prim1(2))*prim1(3)*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
         !dvb(i)=-2*(uspace(i)-prim1(2))*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
        !dvb(i) =-2*(uspace(i)-u)*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
        !dvb(i)=dvh(i)/prim1(3)
     end do
     !do i = 1, unum
      ! print *, "uspace(", i, ") = ", uspace(i)
      ! print *, "weight(", i, ") = ", weight(i)
    ! end do
    end subroutine cal_fv
    
    ! ------------------------------------------------------------
    
    subroutine cal_fv1(cell)
	type(cell_center),intent(in):: cell
	real(kind=8):: du
	integer:: i

	if(associated(dvh)) deallocate(dvh)
	allocate(dvh(unum))
	if(associated(dvb)) deallocate(dvb)
	allocate(dvb(unum))
    

	du=lengthU/unum

	dvh(1)=(-3.0d0*cell%h(1)+4.0d0*cell%h(2)-cell%h(3))/2.0d0/du
	dvb(1)=(-3.0d0*cell%b(1)+4.0d0*cell%b(2)-cell%b(3))/2.0d0/du
    dvh(2)=(cell%h(3)-cell%h(1))/2.0d0/du
	dvb(2)=(cell%b(3)-cell%b(1))/2.0d0/du
    dvh(unum-1)=(cell%h(unum)-cell%h(unum-2))/2.0d0/du
	dvb(unum-1)=(cell%b(unum)-cell%b(unum-2))/2.0d0/du
    dvh(unum)=(3.0d0*cell%h(unum)-4.0d0*cell%h(unum-1)+cell%h(unum-2))/2.0d0/du
    dvb(unum)=(3.0d0*cell%b(unum)-4.0d0*cell%b(unum-1)+cell%b(unum-2))/2.0d0/du

	do i=3,unum-2
		dvh(i)=(-cell%h(i+2)+8.0d0*cell%h(i+1)-8.0d0*cell%h(i-1)+cell%h(i-2))/12.0d0/du
		dvb(i)=(-cell%b(i+2)+8.0d0*cell%b(i+1)-8.0d0*cell%b(i-1)+cell%b(i-2))/12.0d0/du
	end do

    end subroutine cal_fv1
 ! ------------------------------------------------------------
subroutine cal_fv11(cell)
	type(cell_center),intent(in):: cell
    real(kind=8):: prim1(DIM+2)
    real(kind=8):: w1(DIM+2)
	real(kind=8):: du
    real(kind=8), allocatable :: factor(:)
    real(kind=8), allocatable :: hp(:,:)
    real(kind=8), allocatable :: F_term1(:,:)
    real(kind=8), allocatable :: F_term2(:,:)
    real(kind=8), allocatable :: F1(:,:)
    real(kind=8), allocatable :: F2(:,:)
    real(kind=8), allocatable :: v(:)
    real(kind=8), allocatable :: w(:)
    real(kind=8), allocatable :: df(:)
    real(kind=8), allocatable :: dg(:)
    real(kind=8):: rho,u,t
	integer:: i
    integer :: j
    integer :: N
    N=20
	if(associated(dvh)) deallocate(dvh)
	allocate(dvh(unum))
	if(associated(dvb)) deallocate(dvb)
	allocate(dvb(unum))
    if(associated(hH)) deallocate(hH)
	allocate(hH(unum))
	if(associated(bH)) deallocate(bH)
	allocate(bH(unum))
    
    allocate( factor(0:N))
    allocate( df(0:N))
    allocate( dg(0:N))
    allocate( hp(0:N+1,unum))
    allocate( F_term1(unum,0:N))
    allocate( F_term2(unum,0:N))
    allocate( F1(unum,0:N))
    allocate( F2(unum,0:N))
    allocate(v(unum))
    allocate(w(unum))
    
    factor(0) = 1.0
    do i=1,N
        factor(i) = 1.d0
        do j = 1 , i
            factor(i) = factor(i) * j
        end do
    end do
    
    w1(1) = sum(weight*cell%h)
    w1(2) = sum(weight*uspace*cell%h)
    w1(3) = 0.5*(sum(weight*uspace**2*cell%h)+sum(weight*cell%b))
    prim1(1)= w1(1)
    prim1(2)= w1(2)/w1(1)
    prim1(3)= 0.5*w1(1)/(gamma-1.0)/(w1(3)-0.5*w1(2)**2/w1(1))
    
    rho = prim1(1)
    u = prim1(2)
    t = 1.0d0 / prim1(3)
    v = (uspace-u)/dsqrt(t)*dsqrt(2.0d0)
    w = dexp(-v**2/2.d0)/dsqrt(2.d0*PI)
    
    error1=1.0d0
    error2=1.0d0
    do N_order=0,N
     if (error1>0.0005 .or. error2>0.0005) then
       hp(0,:) = 1.d0
       if ( N_order > -1 ) then
         hp(1,:) = v
       end if
      if ( N_order > 0 ) then
        hp(2,:) = v**2 - 1.d0
      end if
     if ( N_order > 1 ) then
       do j=3,N_order+1
         hp(j,:)=v * hp(j-1,:) - ( j - 1.d0 ) * hp(j-2,:)
       end do
     end if
     !print *, "Error1: ", weight
     do j=0,N_order
        df(j)=sum(cell%h * weight * hp(j,:))*dsqrt(2.d0/t)
        dg(j)=sum(cell%b * weight * hp(j,:))*dsqrt(2.d0/t)
         !print *, "df(", j, ") = ", df(j)
         !print *, "dg(", j, ") = ", dg(j)
     end do
    
     F_term1(:,:) = 0.d0
     F_term2(:,:) = 0.d0
     F1(:,:) = 0.d0
     F2(:,:) = 0.d0
     do i=1,unum
        do j=0,N_order
            if (j < 1) then
                 F_term1(i,j) = F_term1(i,j)+df(j)*hp(j+1,i)/factor(j)
                 F_term2(i,j) = F_term2(i,j)+dg(j)*hp(j+1,i)/factor(j)
                 F1(i,j) = F1(i,j)+df(j)*hp(j,i)/factor(j)
                 F2(i,j) = F2(i,j)+dg(j)*hp(j,i)/factor(j)
            end if
            if (j > 0) then
                 F_term1(i,j) = F_term1(i,j-1)+df(j)*hp(j+1,i)/factor(j)
                 F_term2(i,j) = F_term2(i,j-1)+dg(j)*hp(j+1,i)/factor(j)
                 F1(i,j) = F1(i,j-1)+df(j)*hp(j,i)/factor(j)
                 F2(i,j) = F2(i,j-1)+dg(j)*hp(j,i)/factor(j)
            end if
            ! print *, "F_term1(", i, ",", j, ") = ", F_term1(i,j)
        end do
        dvh(i) =-dsqrt(2.0d0/t) * w(i) * F_term1(i,N_order)
        dvb(i) =-dsqrt(2.0d0/t) * w(i) * F_term2(i,N_order)
         !dvh(i)=-2*(uspace(i)-prim1(2))*prim1(3)*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
         !dvb(i)=-2*(uspace(i)-prim1(2))*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
        !dvb(i) =-2*(uspace(i)-u)*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
        !dvb(i)=dvh(i)/prim1(3)
        hH(i) =w(i) * F1(i,N_order)
        bH(i) =w(i) * F2(i,N_order)
        !dvb(i) =-2*(uspace(i)-u)*prim1(1)*(prim1(3)/PI)**(1.0/2.0)*exp(-prim1(3)*(uspace(i)-prim1(2))**2)
        !dvb(i)=dvh(i)/prim1(3)
     end do
     error1=dsqrt(sum(hH-cell%h)**2/sum(cell%h**2))
     error2=dsqrt(sum(bH-cell%b)**2/sum(cell%b**2))
     N_order1=N_order
     end if
      !N_order1=N_order
    end do
    
    
   !N_order1=N_order
   ! print *, "Error1: ", N_order1
    ! print *, "Error2: ", error2
end subroutine cal_fv11
! ------------------------------------------------------------    

 
 

end module GSOURCE
