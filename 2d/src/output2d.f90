module OUTPUT
use ALLDATA
use GRIDMAKER
use GASTHEORY
implicit none
contains
! ------------------------------------------------------------
	subroutine write_macro(loop)
	integer,intent(in):: loop
	character(len=11):: ctemp
	real(kind=8),dimension(:,:,:),allocatable:: csolution
	real(kind=8):: prim(DIM+2)
	integer:: i,j

	allocate(csolution(11,ixmax,iymax))

	do i=1,ixmax
	do j=1,iymax
		prim=get_primary(ctr(i,j)%w)

		csolution(1:3,i,j)=prim(1:3)
		csolution(4,i,j)=1.d0/prim(4)
		csolution(5,i,j)=0.5*prim(1)/prim(4)
		csolution(6:7,i,j)=get_heat_flux(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim)
        csolution(8,i,j)=ctr(i,j)%r
       ! call new_expansion_force_term(ctr(i,j)%h,ctr(i,j)%b, 5,prim, fx, fy, ft1, ft2,error,errort)
        csolution(11,i,j)=1!dsqrt(sum((ctr(i,j)%h-(prim(1)*(prim(4)/PI)*exp(-prim(4)*((uspace*vface(i,j)%cosa+vspace*vface(i,j)%sina-prim(2))**2+(vspace*vface(i,j)%cosa-uspace*vface(i,j)%sina-prim(3))**2))))**2)/sum(ctr(i,j)%h**2))
    end do
    end do
    
    do i=1,ixmax
        do j=1,iymax
        prim=get_primary(ctr(i,1)%w)
        csolution(9,i,j)=ctr(i,1)%r
        csolution(10,i,j)=prim(1)
        end do 
    end do
    

	write(ctemp,'(i11)') loop
    open(unit=fileOutId,file=RSTFILE//trim(adjustl(cTemp))//'.dat',&
         status="replace",action="write")
	write(fileOutId,*) "VARIABLES = X, Y, RHO, U, V, T, P, QX, QY, radius, radius-0, RHO-0,error"
 
	write(fileOutId,'(A,I5,A,I5,A)') 'ZONE I =',grid%xindex, &
        ', J =',grid%yindex, ', DATAPACKING=BLOCK,VARLOCATION=([3-13]=CELLCENTERED)'

	write(fileOutId,"(11(ES23.16,2X))") coordS(1,:,:)
	write(fileOutId,"(11(ES23.16,2X))") coordS(2,:,:)

	do i=1,11
		write(fileOutId,"(11(ES23.16,2X))") csolution(i,:,:)
	end do

	close(fileOutId)
	end subroutine write_macro

! ------------------------------------------------------------
end module OUTPUT