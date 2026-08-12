program test_symmetric_hermite
use CONFIG
use ALLDATA
use VELOCITY
use GASTHEORY
use GSOURCE
use HERMITESYMMETRIC2D
implicit none

integer :: Nm_test
real(8) :: prim(DIM+2), gx, gy, rt, inv_sqrt_rt
real(8) :: err_h, err_b, legacy_err_h, legacy_err_b
real(8) :: symmetric_b_err_h, symmetric_b_err_b
real(8) :: rel_h, rel_b
real(8), allocatable :: h(:,:), b(:,:)
real(8), allocatable :: source_h(:,:), source_b(:,:)
real(8), allocatable :: symmetric_b_h(:,:), symmetric_b_b(:,:)
real(8), allocatable :: legacy_h(:,:), legacy_b(:,:)
real(8), allocatable :: expected_h(:,:), expected_b(:,:)
real(8), allocatable :: ux(:,:), uy(:,:), omega_u(:,:)
real(8), allocatable :: hx2(:,:), hx3(:,:), hx4(:,:)
real(8), allocatable :: hy2(:,:), hy3(:,:), hy4(:,:)

quad_type = 'gauss_hermite'
Tref = 1.d0
nv = 20
call init_velocity()

allocate(h(unum,vnum), b(unum,vnum))
allocate(source_h(unum,vnum), source_b(unum,vnum))
allocate(symmetric_b_h(unum,vnum), symmetric_b_b(unum,vnum))
allocate(legacy_h(unum,vnum), legacy_b(unum,vnum))
allocate(expected_h(unum,vnum), expected_b(unum,vnum))
allocate(ux(unum,vnum), uy(unum,vnum), omega_u(unum,vnum))
allocate(hx2(unum,vnum), hx3(unum,vnum), hx4(unum,vnum))
allocate(hy2(unum,vnum), hy3(unum,vnum), hy4(unum,vnum))

! A manufactured rank-0 through rank-3 symmetric Hermite expansion.
! With T=1 and Rc=0.5, omega(u(v)) is the physicists' GH weight exp(-v^2).
prim = [1.d0, 0.d0, 0.d0, 1.d0]
rt = Rc/prim(4)
inv_sqrt_rt = 1.d0/dsqrt(rt)
ux = (uspace-prim(2))/dsqrt(rt)
uy = (vspace-prim(3))/dsqrt(rt)
omega_u = dexp(-0.5d0*(ux**2+uy**2))/(2.d0*PI)
hx2 = ux**2-1.d0
hx3 = ux**3-3.d0*ux
hx4 = ux**4-6.d0*ux**2+3.d0
hy2 = uy**2-1.d0
hy3 = uy**3-3.d0*uy
hy4 = uy**4-6.d0*uy**2+3.d0

h = omega_u*(1.1d0 + 0.14d0*hx2/2.d0 + 0.07d0*ux*uy + &
    0.03d0*hy3/6.d0)
b = omega_u*(0.8d0 - 0.05d0*ux + 0.09d0*hy2/2.d0 + &
    0.04d0*hx3/6.d0)

gx = 0.13d0
gy = -0.09d0
Nm_test = 3
call symmetric_hermite_force_term_2d(h,b,Nm_test,prim,gx,gy, &
    source_h,source_b,err_h,err_b)
call symmetric_b_hermite_force_term_2d(h,b,Nm_test,prim,gx,gy, &
    symmetric_b_h,symmetric_b_b,symmetric_b_err_h,symmetric_b_err_b)

expected_h = -inv_sqrt_rt*omega_u*( &
    1.1d0*(gx*ux+gy*uy) + &
    0.14d0/2.d0*(gx*hx3+gy*hx2*uy) + &
    0.07d0*(gx*hx2*uy+gy*ux*hy2) + &
    0.03d0/6.d0*(gx*ux*hy3+gy*hy4))
expected_b = -inv_sqrt_rt*omega_u*( &
    0.8d0*(gx*ux+gy*uy) - &
    0.05d0*(gx*hx2+gy*ux*uy) + &
    0.09d0/2.d0*(gx*ux*hy2+gy*hy3) + &
    0.04d0/6.d0*(gx*hx4+gy*hx3*uy))

call require(maxval(abs(source_h-expected_h)) < 1.d-9, &
    'manufactured h force does not match the corrected tensor formula')
call require(maxval(abs(source_b-expected_b)) < 1.d-9, &
    'manufactured b force does not match the corrected tensor formula')
call require(maxval(abs(symmetric_b_h-source_h)) < 1.d-12, &
    'symmetric-b manufactured h force differs from direct-d force')
call require(maxval(abs(symmetric_b_b-source_b)) < 1.d-12, &
    'symmetric-b manufactured b force differs from direct-d force')
call require(err_h < 1.d-9 .and. err_b < 1.d-9, &
    'manufactured distribution is not reconstructed through order 3')
call require(symmetric_b_err_h < 1.d-9 .and. symmetric_b_err_b < 1.d-9, &
    'symmetric-b distribution is not reconstructed through order 3')

! The independent implementation must agree with the existing product-basis
! implementation for a shifted Maxwellian and preserve the force moments.
prim = [1.2d0, 0.23d0, -0.17d0, 1.d0]
call reduced_maxwell(h,b,uspace,vspace,prim)
Nm_test = 6
call check_shifted_maxwell(0.1d0,0.d0)
call check_shifted_maxwell(-0.08d0,0.06d0)

write(*,*) 'PASS: direct-d and normalized symmetric-b Hermite force operators'

contains

subroutine check_shifted_maxwell(ax,ay)
  real(8), intent(in) :: ax, ay
  real(8) :: s0, sx, sy, se

  call symmetric_hermite_force_term_2d(h,b,Nm_test,prim,ax,ay, &
      source_h,source_b,err_h,err_b)
  call symmetric_b_hermite_force_term_2d(h,b,Nm_test,prim,ax,ay, &
      symmetric_b_h,symmetric_b_b,symmetric_b_err_h,symmetric_b_err_b)
  call new_expansion_force_term(h,b,Nm_test,prim,ax,ay, &
      legacy_h,legacy_b,legacy_err_h,legacy_err_b)

  call require(maxval(abs(symmetric_b_h-source_h)) < 1.d-12, &
      'symmetric-b shifted-Maxwellian h force differs from direct-d force')
  call require(maxval(abs(symmetric_b_b-source_b)) < 1.d-12, &
      'symmetric-b shifted-Maxwellian b force differs from direct-d force')
  rel_h = dsqrt(sum((source_h-legacy_h)**2))/max(dsqrt(sum(legacy_h**2)),tiny(1.d0))
  rel_b = dsqrt(sum((source_b-legacy_b)**2))/max(dsqrt(sum(legacy_b**2)),tiny(1.d0))
  call require(rel_h < 1.d-12 .and. rel_b < 1.d-12, &
      'new and legacy Hermite force arrays differ')

  s0 = sum(weight*source_h)
  sx = sum(weight*uspace*source_h) + prim(1)*ax
  sy = sum(weight*vspace*source_h) + prim(1)*ay
  se = 0.5d0*sum(weight*((uspace**2+vspace**2)*source_h+source_b)) &
       + prim(1)*(prim(2)*ax+prim(3)*ay)
  call require(abs(s0) < 2.d-9, 'symmetric Hermite mass moment is inconsistent')
  call require(abs(sx) < 2.d-8 .and. abs(sy) < 2.d-8, &
      'symmetric Hermite momentum moment is inconsistent')
  call require(abs(se) < 5.d-8, 'symmetric Hermite energy moment is inconsistent')
end subroutine check_shifted_maxwell

subroutine require(condition,message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_symmetric_hermite
