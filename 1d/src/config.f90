! ------------------------------------------------------------
! module CONFIG
! Central run-time configuration for the 1D UGKS + Hermite solver.
! All user-tunable parameters live here and are read from a single
! namelist input file (default: input.namelist). Other modules USE
! CONFIG to obtain these values; no parameter is hard-coded elsewhere.
! ------------------------------------------------------------
module CONFIG
implicit none

  ! ---- physical space ----
  integer  :: nx = 600           ! number of physical cells (ixmax)
  real(8)  :: xlength = 1.0d0    ! domain length

  ! ---- velocity space ----
  integer  :: nv = 100           ! number of velocity collocation points (unum)
  real(8)  :: vmax = 5.0d0       ! velocity half-width; domain = [-vmax*sqrt(Tref), +vmax*sqrt(Tref)]
  real(8)  :: Tref = 1.0d0       ! reference temperature for velocity-space scaling
  character(len=20) :: quad_type = 'newton_cotes'   ! 'newton_cotes' or 'gauss_hermite'

  ! ---- body-force term ----
  !   'finite_difference' : velocity-space central difference of f (cal_fv1)
  !   'hermite'            : Hermite spectral expansion of a·∇ᵥf (cal_fv)
  character(len=20) :: force_type = 'hermite'
  integer  :: Nm = 6           ! Hermite expansion order (only used when force_type='hermite')

  ! ---- time stepping ----
  real(8)  :: cfl = 0.5d0
  real(8)  :: maxtime = 0.2d0
  real(8)  :: steady_tol = 1.d-9  ! normalized macro residual (Fourier)
  real(8)  :: mass_tol = 1.d-10   ! relative total-mass drift (Fourier)
  integer  :: steady_steps = 1000 ! consecutive converged steps (Fourier)

  ! ---- case / physics ----
  !   case_id = 2 : 1D hydrostatic equilibrium (well-balanced)
  !   case_id = 3 : Sod shock tube under external force
  !   case_id = 4 : Fourier flow between isothermal plates (diffuse wall)
  integer  :: case_id = 3
  real(8)  :: knudsen = 0.01d0   ! reference Knudsen number
  real(8)  :: phi = -1.0d0       ! external force acceleration
  real(8)  :: TwL = 1.0d0       ! left  wall temperature (diffuse BC, case 4)
  real(8)  :: TwR = 1.1d0       ! right wall temperature (diffuse BC, case 4)
  character(len=40) :: refpfile = ''   ! case 2 reference pressure file (optional)

  ! ---- derived velocity-space bounds (set in read_config) ----
  real(8)  :: umin, umax, du, lengthU

  namelist /params/ nx, xlength, nv, vmax, Tref, quad_type, force_type, Nm, &
                    cfl, maxtime, steady_tol, mass_tol, steady_steps, case_id, &
                    knudsen, phi, TwL, TwR, refpfile

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

    ! derived velocity-space bounds: standard [-vmax, vmax] scaled by sqrt(Tref)
    umin    = -vmax*dsqrt(Tref)
    umax    =  vmax*dsqrt(Tref)
    du      = (umax-umin)/(nv-1)
    lengthU = umax - umin

    call print_config()
  end subroutine read_config
! ------------------------------------------------------------
  subroutine print_config()
    write(*,*) "---- configuration ----"
    write(*,*) "case_id      =", case_id
    write(*,*) "nx           =", nx, "  xlength =", xlength
    write(*,*) "nv           =", nv, "  vmax =", vmax, "  Tref =", Tref
    write(*,*) "quad_type    =", trim(quad_type)
    write(*,*) "force_type   =", trim(force_type), "  Nm =", Nm
    write(*,*) "cfl          =", cfl, "  maxtime =", maxtime
    write(*,*) "steady/mass tol =", steady_tol, mass_tol
    write(*,*) "steady steps =", steady_steps
    write(*,*) "knudsen      =", knudsen, "  phi =", phi
    write(*,*) "TwL/TwR      =", TwL, TwR
    write(*,*) "u-domain     =[", umin, ",", umax, "]  du =", du
    write(*,*) "----------------------"
  end subroutine print_config
! ------------------------------------------------------------
end module CONFIG
