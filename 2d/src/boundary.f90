! ------------------------------------------------------------
! module BOUNDARY
! Periodic boundary condition for the 2D UGKS solver (ghost-cell copy).
! Used by the Poiseuille case (case_id = 5) in the flow (x) direction.
! The specular and diffuse-wall boundary fluxes live in module FLUX
! (calc_flux_specular_vertical/horizon, calc_flux_boundary); this module
! only provides the periodic ghost-cell fill, which has no flux form.
! ------------------------------------------------------------
module BOUNDARY
use ALLDATA
use CONFIG
implicit none

contains
! ------------------------------------------------------------
  subroutine bc_periodic_2d(dirc)
    ! Copy the opposite-end cell into the ghost cell so the internal
    ! flux loop (which spans ixmin..ixmax+1 / iymin..iymax+1) sees a
    ! periodic continuation. No mirroring (unlike specular): a straight
    ! copy of w/h/b/sh/sb.
    integer, intent(in) :: dirc
    integer :: i, j, k, l
    select case(dirc)
    case(Idirc)  ! x direction periodic
      do j = iymin, iymax
        ! left ghost <- rightmost physical cell
        ctr(ixmin-1,j)%w  = ctr(ixmax,j)%w
        ctr(ixmin-1,j)%h  = ctr(ixmax,j)%h
        ctr(ixmin-1,j)%b  = ctr(ixmax,j)%b
        ctr(ixmin-1,j)%sh = ctr(ixmax,j)%sh
        ctr(ixmin-1,j)%sb = ctr(ixmax,j)%sb
        ! right ghost <- leftmost physical cell
        ctr(ixmax+1,j)%w  = ctr(ixmin,j)%w
        ctr(ixmax+1,j)%h  = ctr(ixmin,j)%h
        ctr(ixmax+1,j)%b  = ctr(ixmin,j)%b
        ctr(ixmax+1,j)%sh = ctr(ixmin,j)%sh
        ctr(ixmax+1,j)%sb = ctr(ixmin,j)%sb
      end do
    case(Jdirc)  ! y direction periodic (not used by Poiseuille, kept for symmetry)
      do i = ixmin, ixmax
        ctr(i,iymin-1)%w  = ctr(i,iymax)%w
        ctr(i,iymin-1)%h  = ctr(i,iymax)%h
        ctr(i,iymin-1)%b  = ctr(i,iymax)%b
        ctr(i,iymin-1)%sh = ctr(i,iymax)%sh
        ctr(i,iymin-1)%sb = ctr(i,iymax)%sb
        ctr(i,iymax+1)%w  = ctr(i,iymin)%w
        ctr(i,iymax+1)%h  = ctr(i,iymin)%h
        ctr(i,iymax+1)%b  = ctr(i,iymin)%b
        ctr(i,iymax+1)%sh = ctr(i,iymin)%sh
        ctr(i,iymax+1)%sb = ctr(i,iymin)%sb
      end do
    end select
  end subroutine bc_periodic_2d
! ------------------------------------------------------------
end module BOUNDARY
