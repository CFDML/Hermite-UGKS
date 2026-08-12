module ALLDATA
use CONFIG     ! all run-time parameters (nx, xlength, nv, vmax, Tref, quad_type,
               ! cfl, maxtime, case_id, knudsen, phi, TwL, TwR, refpfile,
               ! umin, umax, du, lengthU) come from here
implicit none

! Math
integer,parameter :: DIM = 1
real(kind=8),parameter :: PI = 3.1415926535d0
real(kind=8),parameter :: smv = tiny(real(1.0,8))
real(kind=8),parameter :: eps = 1.d-6
real(kind=8),parameter :: up = 1.d0
integer,parameter :: RN = 1 !no frame rotation
integer,parameter :: RY = -1 !with frame rotation
integer,parameter :: firstOrder = 1
integer,parameter :: secondOrder = 2
integer,parameter :: thirdOrder = 3
integer,parameter :: fourthOrder = 4
integer,parameter :: centervalue = 1
integer,parameter :: pointsvalue = 2

! Computation (loop state only; cfl/maxtime come from CONFIG)
integer :: interp_method
integer :: iter
integer :: itemp
real(kind=8) :: simtime
real(kind=8) :: dt
real(kind=8) :: res(DIM+2)

! Boundary type (derived from case_id in init_loopSetup, not user-set)
!   bc_type = 1 : specular reflection (cases 2, 3)
!   bc_type = 2 : diffuse / isothermal Maxwellian wall (case 4)
integer :: bc_type

! File
integer,parameter :: fileInId = 17
integer,parameter :: fileOutId = 18
integer,parameter :: fileInitId = 19
integer,parameter :: fileHstId = 20
character(len=40) :: gridfile,inifile,tracefile
character(len=40) :: filename
character(len=13) :: hstfile="shocktube.hst"
character(len=5) :: rstfile="macro"

! Init
integer :: iniType
logical :: isFound,isRestart

! Physics
real(kind=8),parameter :: R0 = 8.314d0
real(kind=8),parameter :: bz = 1.3806488d-23

! Gas property
real(kind=8),parameter :: innerK = 2.d0
real(kind=8),parameter :: Rc = 0.5d0  ! lambda=0.5/T/Rc
real(kind=8),parameter :: NAmass = 39.95d0
real(kind=8),parameter :: molediameter = 4.17d-10
real(kind=8),parameter :: molemass = 6.63d-26
real(kind=8) :: gasR
real(kind=8) :: gamma

! Non-dimensional parameters
real(kind=8),parameter :: prandtl = 0.666666666667d0

! Element
integer :: ixmin,ixmax

! Inifity
real(kind=8) :: rhoInf
real(kind=8) :: velInf
real(kind=8) :: momInf
real(kind=8) :: sosInf
real(kind=8) :: TInf
real(kind=8) :: lambdaInf
real(kind=8) :: EInf
real(kind=8),pointer :: macroVar(:,:) ! Conservative variables

! Reference state
real(kind=8) :: mfp
real(kind=8) :: omega !VHS model in tau
real(kind=8) :: alpha_ref=1.0 !coefficient in HS model
real(kind=8) :: omega_ref=0.5 !coefficient in HS model
real(kind=8) :: muref

! Shock tube
real(kind=8) :: rhoLeft,rhoRight
real(kind=8) :: TLeft,TRight
real(kind=8) :: lambdaLeft,lambdaRight
real(kind=8) :: ELeft,ERight

! Gravity potential field (defined at node points)
real(kind=8),allocatable,dimension(:):: pot

! Velocity space (umin, umax, du, lengthU come from CONFIG; unum set by
! init_velocity in module VELOCITY; uspace/weight allocated there too)
integer :: unum
real(kind=8),allocatable,dimension(:):: uspace
real(kind=8),allocatable,dimension(:):: weight

! Self defined
type:: cell_center
	real(kind=8):: x
	real(kind=8):: length
	real(kind=8):: w(DIM+2)
	real(kind=8):: prim(DIM+2)
	real(kind=8):: sw(DIM+2)
	real(kind=8),allocatable,dimension(:):: h,b
	real(kind=8),allocatable,dimension(:):: sh,sb
end type cell_center

type:: cell_interface
	integer:: bcType
	real(kind=8):: flux(DIM+2)
	real(kind=8),allocatable,dimension(:):: flux_h,flux_b
	real(kind=8):: qf
end type cell_interface

type(cell_center),allocatable,dimension(:):: ctr
type(cell_interface),allocatable,dimension(:):: face

end module ALLDATA
