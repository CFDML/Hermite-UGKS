module OUTPUT
use CONFIG
use ALLDATA
use GRIDMAKER
use GASTHEORY
use GSOURCE
implicit none

    real(kind=8), dimension(:), allocatable :: p_data
    logical :: has_refdata

    contains
! ------------------------------------------------------------
    ! Read reference pressure profile for the hydrostatic case (case 2).
    ! The file origin.dat is an external reference solution; if it is not
    ! present the analytic exp(-x) profile is used instead.
    subroutine read_and_print_p_data(file_name)
    implicit none
    character(len=*), intent(in) :: file_name
    logical :: found
    integer :: fid, i, j, istat, data_index

    has_refdata = .false.
    inquire(file=file_name, exist=found)
    if (.not. found) then
        write(*,*) "reference file '",trim(file_name),"' not found; using exp(-x)"
        return
    end if

    if (allocated(p_data)) deallocate(p_data)
    allocate(p_data(1:600))

    fid = 21
    open(unit=fid, file=file_name, status="old", action="read")
    ! skip header lines
    do i = 1, 270
        read(fid, *)
    end do

    data_index = 0
    do i = 1, 66
        read(fid, *, iostat=istat) (p_data(data_index + j), j = 1, 9)
        if (istat /= 0) then
            write(*,*) 'Error reading reference data or end of file reached.'
            close(fid)
            return
        end if
        data_index = data_index + 9
    end do
    read(fid, *, iostat=istat) (p_data(594 + j), j = 1, 6)
    close(fid)
    has_refdata = .true.

    end subroutine read_and_print_p_data
! ------------------------------------------------------------
    subroutine write_macro(loop)
    integer,intent(in):: loop
    character(len=8):: ctemp
    real(kind=8),dimension(:,:),allocatable:: csolution
    real(kind=8):: prim(DIM+2)
    real(kind=8):: pref
    integer :: i

    allocate(csolution(9,ixmax))

    select case(case_id)
    case(2)
        call read_and_print_p_data(refpfile)
    end select

    do i=ixmin,ixmax
        prim=get_primary(ctr(i)%w)
        csolution(1:2,i)=prim(1:2)
        csolution(3,i)=get_temperature_from_lambda(prim(3))
        csolution(4,i)=Rc*csolution(1,i)*csolution(3,i)
        ! reference pressure and perturbation (P - P_ref)
        if (case_id==2 .and. has_refdata) then
            pref = p_data(i)
        else
            pref = exp(-ctr(i)%x)
        end if
        csolution(5,i)=csolution(4,i)-pref
        csolution(6,i)=get_heat_flux(ctr(i)%h,ctr(i)%b,prim)
        csolution(7,i)=error1
        csolution(8,i)=error2
        csolution(9,i)=N_order1
    end do

    write(ctemp,'(i8)') loop
    open(unit=fileOutId,file=RSTFILE//trim(adjustl(cTemp))//'.dat',&
         status="replace",action="write")
    write(fileOutId,*) "VARIABLES = X, RHO, U, T, P, PERTURBATION, Q, error1, error2, N"

    write(fileOutId,'(A,I5,A,I5,A)') 'ZONE I =',ixmax-ixmin+1,', DATAPACKING=BLOCK'
    write(fileOutId,"(9(ES23.16,2X))") ctr(ixmin:ixmax)%x
    do i=1,9
        write(fileOutId,"(9(ES23.16,2X))") csolution(i,:)
    end do
    close(fileOutId)

    if (allocated(p_data)) deallocate(p_data)
    deallocate(csolution)
    end subroutine write_macro
! ------------------------------------------------------------
end module OUTPUT
