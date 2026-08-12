module ALLDATA
use CONFIG     ! all run-time parameters (nx,ny,xlength,ylength,nv,vmax,Tref,
               ! quad_type,cfl,maxtime,case_id,knudsen,phi,TwL,TwR,umin,umax,du,
               ! vmin,vmax_,dv,lengthU,lengthV) come from here
implicit none

! Math
integer,parameter :: DIM = 2
real(kind=8),parameter :: PI = 3.1415926535d0
real(kind=8),parameter :: smv = 1.d-8
real(kind=8),parameter :: eps = 1.d-7
real(kind=8),parameter :: up = 1.d0
integer,parameter :: Idirc = 1
integer,parameter :: Jdirc = 2
integer,parameter :: RN = 1 !no frame rotation
integer,parameter :: RY = -1 !with frame rotation
integer,parameter :: firstOrder = 1
integer,parameter :: secondOrder = 2
integer,parameter :: thirdOrder = 3
integer,parameter :: fourthOrder = 4
integer,parameter :: centervalue = 1
integer,parameter :: pointsvalue = 2

! Computation (loop state; cfl/maxtime come from CONFIG)
integer :: interp_method
integer :: iter
integer :: itemp
real(kind=8) :: simtime
real(kind=8) :: time1,time2,time3,time4
real(kind=8) :: dt
real(kind=8) :: res(DIM+2)

! File
integer,parameter :: fileInId = 17
integer,parameter :: fileOutId = 18
integer,parameter :: fileInitId = 19
integer,parameter :: fileHstId = 20
character(len=40) :: gridfile,inifile,tracefile
character(len=40) :: filename
character(len=13) :: hstfile="rtcompute.hst"
character(len=5) :: rstfile="macro"
character(len=5) :: mstfile="micro"

! Init
integer :: iniType
logical :: isFound,isRestart

! Physics
real(kind=8),parameter :: R0 = 8.314d0
real(kind=8),parameter :: bz = 1.3806488d-23

! Gas property (Ar gas)
real(kind=8),parameter :: innerK = 1.d0
real(kind=8),parameter :: Rc = 0.5d0  ! lambda=0.5/T/Rc
real(kind=8),parameter :: NAmass = 39.95d0
real(kind=8),parameter :: molediameter = 4.17d-10
real(kind=8),parameter :: molemass = 6.63d-26
real(kind=8) :: gasR
real(kind=8) :: gamma

! Non-dimensional parameters
! Case dependent: Pr=1 for Poiseuille and Pr=2/3 for Rayleigh--Taylor.
real(kind=8) :: prandtl = 0.666666666667d0
real(kind=8),parameter :: mach = 0.15d0

! Element
integer :: npe !node per element
integer :: ixmin,ixmax
integer :: iymin,iymax

! Inifity
real(kind=8) :: velInf(2)
real(kind=8) :: momInf(2)
real(kind=8) :: sosInf,lambdaInf
real(kind=8) :: rhoInf
real(kind=8) :: TInf
real(kind=8) :: EInf

! Reference state
real(kind=8) :: mfp
real(kind=8) :: omega !VHS model in tau
real(kind=8) :: alpha_ref=1.0 !coefficient in HS model
real(kind=8) :: omega_ref=0.5 !coefficient in HS model
real(kind=8) :: muref

! Gravity convection (phi comes from CONFIG)
real(kind=8):: radius
real(kind=8):: bc(4)=(/1.0,0.0,0.0,1.6666667/)

! velocity space (umin,umax,du,vmin,vmax_,dv,lengthU,lengthV from CONFIG;
! unum/vnum set by init_velocity; uspace/vspace/weight allocated there)
integer:: unum,vnum
real(kind=8),allocatable,dimension(:,:):: uspace,vspace
real(kind=8),allocatable,dimension(:,:):: weight

! self defined
type:: cell_center
	real(kind=8):: x,y
	real(kind=8):: r,theta
	real(kind=8):: area
	real(kind=8):: length(DIM)
	real(kind=8):: w(DIM+2)
	real(kind=8):: prim(DIM+2)
	real(kind=8):: sw(DIM+2,2)
	real(kind=8),allocatable,dimension(:,:):: h,b
	real(kind=8),allocatable,dimension(:,:,:):: sh,sb
end type cell_center

type:: cell_interface
	integer:: bcType
	real(kind=8):: length
	real(kind=8):: cosa,sina
	real(kind=8):: theta
	real(kind=8):: flux(DIM+2)
	real(kind=8),allocatable,dimension(:,:):: flux_h,flux_b
	real(kind=8):: qf(DIM)
end type cell_interface

type(cell_center),allocatable,dimension(:,:):: ctr
type(cell_interface),allocatable,dimension(:,:):: vface,hface

end module ALLDATA
