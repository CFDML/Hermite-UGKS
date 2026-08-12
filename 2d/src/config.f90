! ------------------------------------------------------------
! module CONFIG
! Central run-time configuration for the 2D UGKS + Hermite solver.
! All user-tunable parameters live here and are read from a single
! namelist input file (default: input.namelist).
! ------------------------------------------------------------
module CONFIG
implicit none

  ! ---- physical space ----
  integer  :: nx = 3            ! number of cells in x
  integer  :: ny = 100          ! number of cells in y
  real(8)  :: xlength = 0.1d0   ! domain length in x
  real(8)  :: ylength = 1.0d0   ! domain length in y

  ! ---- velocity space ----
  integer  :: nv = 16           ! number of velocity points per direction (unum=vnum)
  real(8)  :: vmax = 5.0d0      ! velocity half-width (domain = ±vmax*sqrt(Tref))
  real(8)  :: Tref = 1.0d0      ! reference temperature for velocity scaling
  character(len=20) :: quad_type = 'gauss_hermite'   ! 'newton_cotes' or 'gauss_hermite'

  ! ---- body-force term ----
  !   'finite_difference' : velocity-space central difference of f (cal_fsource1)
  !   'hermite'            : Hermite spectral expansion of a·∇ᵥf (new_expansion_force_term)
  !   'hermite_symmetric'   : corrected direct-d symmetric Hermite expansion
  !   'hermite_symmetric_b' : corrected symmetric-b Hermite expansion
  character(len=20) :: force_type = 'hermite_symmetric'
  integer  :: Nm = 7           ! Hermite order for either Hermite force operator

  ! ---- time stepping ----
  real(8)  :: cfl = 0.8d0
  real(8)  :: maxtime = 200.d0

  ! ---- case / physics ----
  !   case_id = 6 : Rayleigh-Taylor instability (quarter domain)
  !   case_id = 5 : Poiseuille flow (x-periodic, y-diffuse walls)
  integer  :: case_id = 5
  real(8)  :: knudsen = 2.d-2
  real(8)  :: phi = 0.1d0       ! external force (RT: -1.5 gravity; Poiseuille: +0.1 in x)
  real(8)  :: TwL = 1.0d0      ! lower/left  wall temperature (diffuse BC)
  real(8)  :: TwR = 1.0d0      ! upper/right wall temperature (diffuse BC)

  ! ---- derived velocity-space bounds (set in read_config) ----
  real(8)  :: umin, umax, du
  real(8)  :: vmin, vmax_, dv, lengthU, lengthV

  namelist /params/ nx, ny, xlength, ylength, nv, vmax, Tref, quad_type, &
                    force_type, Nm, cfl, maxtime, case_id, knudsen, phi, TwL, TwR

contains
! ------------------------------------------------------------
  subroutine read_config(fname)
    character(len=*), intent(in) :: fname
    integer :: u, ios
    logical :: found

    inquire(file=fname, exist=found)
    if (.not. found) then
      write(*,*) "config file '",trim(fname),"' not found; using built-in defaults"
    else
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) then
        write(*,*) "error opening '",trim(fname),"'; using defaults"
      else
        read(u, nml=params, iostat=ios)
        if (ios /= 0) write(*,*) "warning: namelist parse error in '",trim(fname), &
                                  "'; using partial defaults"
        close(u)
      end if
    end if

    ! derived velocity-space bounds: [-vmax,vmax]*sqrt(Tref), same in u and v
    umin    = -vmax*dsqrt(Tref)
    umax    =  vmax*dsqrt(Tref)
    du      = (umax-umin)/(nv-1)
    lengthU = umax - umin
    vmin    = umin
    vmax_   = umax
    dv      = du
    lengthV = lengthU

    call print_config()
  end subroutine read_config
! ------------------------------------------------------------
  subroutine print_config()
    write(*,*) "---- 2D configuration ----"
    write(*,*) "case_id      =", case_id
    write(*,*) "nx,ny        =", nx, ny, "  xlength,ylength =", xlength, ylength
    write(*,*) "nv           =", nv, "  vmax =", vmax, "  Tref =", Tref
    write(*,*) "quad_type    =", trim(quad_type)
    write(*,*) "force_type   =", trim(force_type), "  Nm =", Nm
    write(*,*) "cfl          =", cfl, "  maxtime =", maxtime
    write(*,*) "knudsen      =", knudsen, "  phi =", phi
    write(*,*) "TwL/TwR      =", TwL, TwR
    write(*,*) "u-domain     =[", umin, ",", umax, "]  du =", du
    write(*,*) "-------------------------"
  end subroutine print_config
! ------------------------------------------------------------
end module CONFIG
