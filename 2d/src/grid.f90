module GRIDMAKER
use CONFIG
use ALLDATA
implicit none

real(kind=8),pointer:: coordUS(:,:) !非结构网格
real(kind=8),pointer:: coordS(:,:,:) !结构网格

type gridType
	integer:: numNp,numElem,nodePerElem,numEdge
	! mapping between Elem and Np
	integer,pointer:: mapEN(:,:) !每个element对应的四个node points
	integer,pointer:: edgeMapNp(:,:) !每条边对应两个节点
	integer,pointer:: edgeNeighborElem(:,:) !每条边对应两个单元
	integer,pointer:: elemNeighbor(:,:) !单元相邻的四个单元
	integer,pointer:: elemMapEdge(:,:)
	! index of structured node points/elements
	integer:: xindex,yindex
	! map from 2D list to 1D list
	integer,pointer:: mapStoUns(:,:)
end type gridType

type(gridType):: grid	

integer,parameter:: selfgnt = 1
integer,parameter:: loadmsh = 2

contains
! ------------------------------------------------------------
	subroutine gridCheck()
	logical:: found

	inquire(file=gridfile,exist=found)
	if (.not. found) then
		write(*,*) 'Grid file "',trim(gridfile),'" does not exist'
		if(initype==loadmsh) stop
	else
		write(*,*) 'Grid file is ',trim(gridfile)
	end if
	end subroutine gridCheck
! ------------------------------------------------------------
	subroutine gridInit()

	select case(initype)
		case(selfgnt)
			write(*,*) "Generating grid"
			call gntSGrid !人为设置网格参数
		case(loadmsh)
			write(*,*) "Loading grid"
	end select
	end subroutine gridInit
! ------------------------------------------------------------
	subroutine gntSGrid
	real(kind=8):: dx,dy
	integer:: i,j

	write(*,*) "The grid is structured"
	! gridType: node count = cell count + 1
	grid%xindex=nx+1
	grid%yindex=ny+1
	grid%numNp=grid%xindex*grid%yindex
	grid%numElem=(grid%xindex-1)*(grid%yindex-1)
	grid%nodePerElem=4

	! element
	ixmin=1
	ixmax=grid%xindex-1
	iymin=1
	iymax=grid%yindex-1

	if(associated(coordS)) then
		deallocate(coordS)
	endif
	allocate(coordS(2,grid%xIndex,grid%yIndex))
	
	dx=xlength/(ixmax-ixmin+1)
	dy=ylength/(iymax-iymin+1)
	do i=1,grid%xindex
	do j=1,grid%yindex
		coordS(1,i,j)=(i-1)*dx
		coordS(2,i,j)=(j-1)*dy
	end do
	end do

	write(*,*) "xlength: ",xlength
	write(*,*) "ylength: ",ylength
	write(*,*) "x cell number: ",ixmax-ixmin+1
	write(*,*) "y cell number: ",iymax-iymin+1

	return
	end subroutine gntSGrid
! ------------------------------------------------------------
end module GRIDMAKER