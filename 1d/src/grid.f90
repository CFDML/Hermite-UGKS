module GRIDMAKER
use CONFIG
use ALLDATA
implicit none

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
	real(kind=8)::  dx,dy
	integer:: i,j

	write(*,*) "The grid is structured"

	! element
	ixmin=1
	ixmax=nx
	
	write(*,*) "xlength: ",xlength
	write(*,*) "x cell number: ",ixmax-ixmin+1
	write(*,*) "ghost cell: ",0
	
	return
	end subroutine gntSGrid
! ------------------------------------------------------------
end module GRIDMAKER