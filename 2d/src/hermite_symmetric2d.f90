! ------------------------------------------------------------
! Correct two-dimensional symmetric-tensor Hermite force operator.
!
! For n=r+s, the independent tensor component d_(r,s) represents all
! permutations with r x-indices and s y-indices.  The factor
!
!   (1/n!) * binomial(n,r) = 1/(r! s!)
!
! converts the full symmetric-tensor contraction to triangular storage.
! Both the direct d^(n) contraction and the equivalent normalized
! b^(n+1)=Sym(a tensor-product d^(n)) contraction are provided.
! ------------------------------------------------------------
module HERMITESYMMETRIC2D
use ALLDATA
implicit none
private

public :: symmetric_hermite_force_term_2d
public :: symmetric_b_hermite_force_term_2d

contains

subroutine symmetric_hermite_force_term_2d(h,b,order,prim,g_x,g_y, &
    force_h,force_b,error_h,error_b)
  integer, intent(in) :: order
  real(8), intent(in) :: h(unum,vnum), b(unum,vnum)
  real(8), intent(in) :: prim(DIM+2), g_x, g_y
  real(8), intent(out) :: force_h(unum,vnum), force_b(unum,vnum)
  real(8), intent(out) :: error_h, error_b

  integer :: n, r, s
  real(8) :: rt, inv_sqrt_rt, jacobian, denominator
  real(8), allocatable :: u_x(:,:), u_y(:,:), omega_u(:,:)
  real(8), allocatable :: hermite_x(:,:,:), hermite_y(:,:,:)
  real(8), allocatable :: factorial(:)
  real(8), allocatable :: coefficient_h(:,:), coefficient_b(:,:)
  real(8), allocatable :: reconstruction_h(:,:), reconstruction_b(:,:)

  if (order < 0) error stop 'symmetric Hermite: order must be non-negative'
  if (prim(4) <= 0.d0) error stop 'symmetric Hermite: temperature must be positive'

  ! prim(4)=lambda=1/T and Rc=R=0.5 in this solver.
  rt = Rc/prim(4)
  inv_sqrt_rt = 1.d0/dsqrt(rt)
  jacobian = 1.d0/rt

  allocate(u_x(unum,vnum),u_y(unum,vnum),omega_u(unum,vnum))
  allocate(hermite_x(0:order+1,unum,vnum))
  allocate(hermite_y(0:order+1,unum,vnum))
  allocate(factorial(0:order))
  allocate(coefficient_h(0:order,0:order))
  allocate(coefficient_b(0:order,0:order))
  allocate(reconstruction_h(unum,vnum),reconstruction_b(unum,vnum))

  u_x = (uspace-prim(2))/dsqrt(rt)
  u_y = (vspace-prim(3))/dsqrt(rt)
  omega_u = dexp(-0.5d0*(u_x**2+u_y**2))/(2.d0*PI)

  call build_probabilists_hermite(u_x,order+1,hermite_x)
  call build_probabilists_hermite(u_y,order+1,hermite_y)

  factorial(0) = 1.d0
  do n = 1, order
    factorial(n) = dble(n)*factorial(n-1)
  end do

  coefficient_h = 0.d0
  coefficient_b = 0.d0
  do n = 0, order
    do r = 0, n
      s = n-r
      coefficient_h(r,s) = jacobian*sum(weight*h*hermite_x(r,:,:)*hermite_y(s,:,:))
      coefficient_b(r,s) = jacobian*sum(weight*b*hermite_x(r,:,:)*hermite_y(s,:,:))
    end do
  end do

  force_h = 0.d0
  force_b = 0.d0
  reconstruction_h = 0.d0
  reconstruction_b = 0.d0
  do n = 0, order
    do r = 0, n
      s = n-r
      denominator = factorial(r)*factorial(s)

      reconstruction_h = reconstruction_h + omega_u*coefficient_h(r,s) * &
          hermite_x(r,:,:)*hermite_y(s,:,:)/denominator
      reconstruction_b = reconstruction_b + omega_u*coefficient_b(r,s) * &
          hermite_x(r,:,:)*hermite_y(s,:,:)/denominator

      force_h = force_h-inv_sqrt_rt*omega_u*coefficient_h(r,s)/denominator * &
          (g_x*hermite_x(r+1,:,:)*hermite_y(s,:,:) + &
           g_y*hermite_x(r,:,:)*hermite_y(s+1,:,:))
      force_b = force_b-inv_sqrt_rt*omega_u*coefficient_b(r,s)/denominator * &
          (g_x*hermite_x(r+1,:,:)*hermite_y(s,:,:) + &
           g_y*hermite_x(r,:,:)*hermite_y(s+1,:,:))
    end do
  end do

  error_h = dsqrt(sum((h-reconstruction_h)**2)) / &
      max(dsqrt(sum(h**2)),tiny(1.d0))
  error_b = dsqrt(sum((b-reconstruction_b)**2)) / &
      max(dsqrt(sum(b**2)),tiny(1.d0))
end subroutine symmetric_hermite_force_term_2d

! ------------------------------------------------------------
! Equivalent force operator written through the independent components
! of the normalized fully symmetric tensor
!
!   b^(n+1) = Sym(a tensor-product d^(n)).
!
! For r+s=n+1,
!
!   b_(r,s) = [r*a_x*d_(r-1,s) + s*a_y*d_(r,s-1)]/(n+1),
!
! where a term with a negative index is omitted.  Combining the full
! symmetric contraction with the 1/n! coefficient gives the independent
! component factor (n+1)/(r!s!).
! ------------------------------------------------------------
subroutine symmetric_b_hermite_force_term_2d(h,b,order,prim,g_x,g_y, &
    force_h,force_b,error_h,error_b)
  integer, intent(in) :: order
  real(8), intent(in) :: h(unum,vnum), b(unum,vnum)
  real(8), intent(in) :: prim(DIM+2), g_x, g_y
  real(8), intent(out) :: force_h(unum,vnum), force_b(unum,vnum)
  real(8), intent(out) :: error_h, error_b

  integer :: n, r, s
  real(8) :: rt, inv_sqrt_rt, jacobian, denominator
  real(8), allocatable :: u_x(:,:), u_y(:,:), omega_u(:,:)
  real(8), allocatable :: hermite_x(:,:,:), hermite_y(:,:,:)
  real(8), allocatable :: factorial(:)
  real(8), allocatable :: d_h(:,:), d_b(:,:)
  real(8), allocatable :: symmetric_b_h(:,:), symmetric_b_b(:,:)
  real(8), allocatable :: reconstruction_h(:,:), reconstruction_b(:,:)

  if (order < 0) error stop 'symmetric-b Hermite: order must be non-negative'
  if (prim(4) <= 0.d0) error stop 'symmetric-b Hermite: temperature must be positive'

  ! prim(4)=lambda=1/T and Rc=R=0.5 in this solver.
  rt = Rc/prim(4)
  inv_sqrt_rt = 1.d0/dsqrt(rt)
  jacobian = 1.d0/rt

  allocate(u_x(unum,vnum),u_y(unum,vnum),omega_u(unum,vnum))
  allocate(hermite_x(0:order+1,unum,vnum))
  allocate(hermite_y(0:order+1,unum,vnum))
  allocate(factorial(0:order+1))
  allocate(d_h(0:order,0:order),d_b(0:order,0:order))
  allocate(symmetric_b_h(0:order+1,0:order+1))
  allocate(symmetric_b_b(0:order+1,0:order+1))
  allocate(reconstruction_h(unum,vnum),reconstruction_b(unum,vnum))

  u_x = (uspace-prim(2))/dsqrt(rt)
  u_y = (vspace-prim(3))/dsqrt(rt)
  omega_u = dexp(-0.5d0*(u_x**2+u_y**2))/(2.d0*PI)

  call build_probabilists_hermite(u_x,order+1,hermite_x)
  call build_probabilists_hermite(u_y,order+1,hermite_y)

  factorial(0) = 1.d0
  do n = 1, order+1
    factorial(n) = dble(n)*factorial(n-1)
  end do

  d_h = 0.d0
  d_b = 0.d0
  reconstruction_h = 0.d0
  reconstruction_b = 0.d0
  do n = 0, order
    do r = 0, n
      s = n-r
      denominator = factorial(r)*factorial(s)
      d_h(r,s) = jacobian*sum(weight*h*hermite_x(r,:,:)*hermite_y(s,:,:))
      d_b(r,s) = jacobian*sum(weight*b*hermite_x(r,:,:)*hermite_y(s,:,:))
      reconstruction_h = reconstruction_h + omega_u*d_h(r,s) * &
          hermite_x(r,:,:)*hermite_y(s,:,:)/denominator
      reconstruction_b = reconstruction_b + omega_u*d_b(r,s) * &
          hermite_x(r,:,:)*hermite_y(s,:,:)/denominator
    end do
  end do

  symmetric_b_h = 0.d0
  symmetric_b_b = 0.d0
  force_h = 0.d0
  force_b = 0.d0
  do n = 0, order
    do r = 0, n+1
      s = n+1-r

      if (r > 0) then
        symmetric_b_h(r,s) = symmetric_b_h(r,s) + &
            dble(r)*g_x*d_h(r-1,s)
        symmetric_b_b(r,s) = symmetric_b_b(r,s) + &
            dble(r)*g_x*d_b(r-1,s)
      end if
      if (s > 0) then
        symmetric_b_h(r,s) = symmetric_b_h(r,s) + &
            dble(s)*g_y*d_h(r,s-1)
        symmetric_b_b(r,s) = symmetric_b_b(r,s) + &
            dble(s)*g_y*d_b(r,s-1)
      end if
      symmetric_b_h(r,s) = symmetric_b_h(r,s)/dble(n+1)
      symmetric_b_b(r,s) = symmetric_b_b(r,s)/dble(n+1)

      denominator = factorial(r)*factorial(s)
      force_h = force_h - inv_sqrt_rt*omega_u*dble(n+1) * &
          symmetric_b_h(r,s)*hermite_x(r,:,:)*hermite_y(s,:,:)/denominator
      force_b = force_b - inv_sqrt_rt*omega_u*dble(n+1) * &
          symmetric_b_b(r,s)*hermite_x(r,:,:)*hermite_y(s,:,:)/denominator
    end do
  end do

  error_h = dsqrt(sum((h-reconstruction_h)**2)) / &
      max(dsqrt(sum(h**2)),tiny(1.d0))
  error_b = dsqrt(sum((b-reconstruction_b)**2)) / &
      max(dsqrt(sum(b**2)),tiny(1.d0))
end subroutine symmetric_b_hermite_force_term_2d

subroutine build_probabilists_hermite(u,max_order,basis)
  integer, intent(in) :: max_order
  real(8), intent(in) :: u(:,:)
  real(8), intent(out) :: basis(0:max_order,size(u,1),size(u,2))
  integer :: n

  basis(0,:,:) = 1.d0
  if (max_order == 0) return

  basis(1,:,:) = u
  do n = 2, max_order
    basis(n,:,:) = u*basis(n-1,:,:)-dble(n-1)*basis(n-2,:,:)
  end do
end subroutine build_probabilists_hermite

end module HERMITESYMMETRIC2D
