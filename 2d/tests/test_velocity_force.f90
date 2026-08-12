program test_velocity_force
use CONFIG
use ALLDATA
use VELOCITY
use GASTHEORY
use GSOURCE
implicit none

integer, parameter :: orders(3) = [16, 18, 20]
integer :: k
real(8) :: mass, mx, my, mxx, myy, tol
real(8) :: prim(DIM+2), err_h, err_b
real(8) :: tau_expected
real(8), allocatable :: h(:,:), b(:,:), source_h(:,:), source_b(:,:)

quad_type = 'gauss_hermite'
Tref = 1.d0
! Nv=16 is the paper setting for Kn=0.02.  Plain-integral GH weights
! approximate the non-weighted Maxwellian moments spectrally; at Nv=16 the
! largest tested second-moment truncation error is about 7.3e-7.
tol = 1.d-6

do k = 1, size(orders)
  nv = orders(k)
  call clear_velocity()
  call init_velocity()

  call require(all(uspace(2:,1) > uspace(:unum-1,1)), &
               'GH nodes are not strictly increasing')
  call require(maxval(abs(uspace(:,1) + uspace(unum:1:-1,1))) < 2.d-13, &
               'GH x nodes are not symmetric')
  call require(maxval(abs(vspace(1,:) + vspace(1,vnum:1:-1))) < 2.d-13, &
               'GH y nodes are not symmetric')
  call require(all(weight > 0.d0), 'GH weights are not positive')

  mass = sum(weight*exp(-(uspace**2+vspace**2)/2.d0)/(2.d0*PI))
  mx = sum(weight*uspace*exp(-(uspace**2+vspace**2)/2.d0)/(2.d0*PI))
  my = sum(weight*vspace*exp(-(uspace**2+vspace**2)/2.d0)/(2.d0*PI))
  mxx = sum(weight*uspace**2*exp(-(uspace**2+vspace**2)/2.d0)/(2.d0*PI))
  myy = sum(weight*vspace**2*exp(-(uspace**2+vspace**2)/2.d0)/(2.d0*PI))
  call require(abs(mass-1.d0) < tol, '2D Maxwell mass is inaccurate')
  call require(abs(mx) < tol .and. abs(my) < tol, '2D Maxwell first moment is inaccurate')
  call require(abs(mxx-1.d0) < tol .and. abs(myy-1.d0) < tol, &
               '2D Maxwell second moment is inaccurate')
end do

prim = [1.2d0, 0.d0, 0.d0, 0.8d0]
case_id = 5
muref = 0.5d0*sqrt(PI)*0.02d0
call require(abs(get_tau(prim)-muref/prim(1)) < 2.d-15, &
             'Poiseuille collision time does not match its reference model')
case_id = 6
omega = 0.72d0
tau_expected = muref*2.d0*prim(4)**(1.d0-omega)/prim(1)
call require(abs(get_tau(prim)-tau_expected) < 2.d-15, &
             'RT collision time does not match its VHS model')

nv = 20
Nm = 6
call clear_velocity()
call init_velocity()
allocate(h(unum,vnum), b(unum,vnum), source_h(unum,vnum), source_b(unum,vnum))
prim = [1.2d0, 0.23d0, -0.17d0, 1.d0]
call reduced_maxwell(h,b,uspace,vspace,prim)

call check_force(0.1d0, 0.d0)
call check_force(-0.08d0, 0.06d0)

write(*,*) 'PASS: 2D Gauss-Hermite and Hermite force moments'

contains

subroutine clear_velocity()
  if (allocated(uspace)) deallocate(uspace)
  if (allocated(vspace)) deallocate(vspace)
  if (allocated(weight)) deallocate(weight)
end subroutine clear_velocity

subroutine check_force(gx,gy)
  real(8), intent(in) :: gx,gy
  real(8) :: s0,sx,sy,se
  call new_expansion_force_term(h,b,Nm,prim,gx,gy,source_h,source_b,err_h,err_b)
  s0 = sum(weight*source_h)
  sx = sum(weight*uspace*source_h) + prim(1)*gx
  sy = sum(weight*vspace*source_h) + prim(1)*gy
  se = 0.5d0*sum(weight*((uspace**2+vspace**2)*source_h+source_b)) &
       + prim(1)*(prim(2)*gx+prim(3)*gy)
  call require(abs(s0) < 2.d-9, 'Hermite force mass moment is inconsistent')
  call require(abs(sx) < 2.d-8 .and. abs(sy) < 2.d-8, &
               'Hermite force momentum moment is inconsistent')
  call require(abs(se) < 5.d-8, 'Hermite force energy moment is inconsistent')
end subroutine check_force

subroutine require(condition,message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message), ' nv=', nv
    error stop 1
  end if
end subroutine require

end program test_velocity_force
