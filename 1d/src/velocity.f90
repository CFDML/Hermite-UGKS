! ------------------------------------------------------------
! module VELOCITY
! Velocity-space node/weight generation for the 1D UGKS solver.
! Two quadratures, selected by quad_type from module CONFIG:
!   'newton_cotes'  : uniform nodes, 6th-order closed Newton-Cotes weights
!   'gauss_hermite' : roots of the physicist's Hermite polynomial H_nv,
!                     weights from Eq.(34) of arXiv:2507.10021
! Both produce (uspace, weight) on the physical velocity axis, i.e. the
! weights are plain quadrature weights for ∫ g(v) dv (no extra e^{v^2}
! rescaling). Nodes are placed on [-vmax*sqrt(Tref), +vmax*sqrt(Tref)]
! for Newton-Cotes; Gauss-Hermite nodes are the standard roots scaled by
! sqrt(Tref) (so the Maxwellian e^{-v^2/(2Tref)} is integrated exactly).
! ------------------------------------------------------------
module VELOCITY
use ALLDATA
use CONFIG
implicit none

contains
! ------------------------------------------------------------
  subroutine init_velocity()
    ! dispatcher: choose quadrature by quad_type
    select case(trim(quad_type))
    case('newton_cotes')
      call init_velocity_newton()
    case('gauss_hermite')
      call init_velocity_gauss_hermite()
    case default
      write(*,*) "unknown quad_type '",trim(quad_type),"'; using newton_cotes"
      call init_velocity_newton()
    end select
  end subroutine init_velocity
! ------------------------------------------------------------
  subroutine init_velocity_newton()
    integer :: i
    ! unum, umin, umax, du come from CONFIG (unum=nv)
    unum = nv
    allocate(uspace(unum))
    allocate(weight(unum))

    du = (umax-umin)/(unum-1)
    do i = 1, unum
      uspace(i) = umin + (i-1)*du
      weight(i) = newton_coeff(i,unum)*du
    end do

    write(*,*) "velocity: Newton-Cotes, unum =", unum, " du =", du
  end subroutine init_velocity_newton
! ------------------------------------------------------------
  pure function newton_coeff(idx,num)
    ! 6th-order closed Newton-Cotes coefficient (multiplies du)
    integer, intent(in) :: idx, num
    real(8) :: newton_coeff
    if (idx==1 .or. idx==num) then
      newton_coeff = 14.0d0/45.0d0
    else if (mod(idx-5,4)==0) then
      newton_coeff = 28.0d0/45.0d0
    else if (mod(idx-3,4)==0) then
      newton_coeff = 24.0d0/45.0d0
    else
      newton_coeff = 64.0d0/45.0d0
    end if
  end function newton_coeff
! ------------------------------------------------------------
  subroutine init_velocity_gauss_hermite()
    ! Gauss-Hermite quadrature: nodes = roots of physicist's H_n(v),
    ! weights w_j = 2^(n-1) n! sqrt(pi) / (n^2 [H'_{n-1}(v_j)]^2).
    ! Physicist's H_n satisfies  ∫ H_m H_n e^{-v^2} dv = sqrt(pi) 2^n n! delta_mn.
    ! The standard GH rule integrates ∫ e^{-v^2} I(v) dv = sum w_j I(v_j).
    ! To integrate a plain function g(v) = ∫ g(v) dv we store the weight as
    ! w_j * exp(v_j^2) so that  sum_j [w_j exp(v_j^2)] g(v_j) = ∫ g(v) dv.
    ! Nodes are then scaled by sqrt(Tref) so the physical Maxwellian is exact.
    integer :: n
    integer :: i, j, iter, m
    real(8) :: z, zold, p0, p1, p2, dp, scale
    real(8), allocatable :: roots(:), wgh(:), positive_roots(:), positive_weights(:)
    real(8), parameter :: eps_tol = 1.d-14

    n = nv
    if (n < 2) then
      write(*,*) "gauss_hermite: nv must be >= 2, got ", n
      stop
    end if

    m = (n + 1)/2
    allocate(roots(n), wgh(n), positive_roots(m), positive_weights(m))
    scale = dsqrt(Tref)

    ! Find the positive roots from the outside inward.  The extrapolated
    ! starting values keep each Newton iteration in the basin of a different
    ! root.  A normalized recurrence avoids overflow in H_n at useful orders.
    z = 0.d0
    do i = 1, m
      if (i == 1) then
        z = dsqrt(2.d0*n + 1.d0) - 1.85575d0*(2.d0*n + 1.d0)**(-1.d0/6.d0)
      else if (i == 2) then
        z = z - 1.14d0*dble(n)**0.426d0/z
      else if (i == 3) then
        z = 1.86d0*z - 0.86d0*positive_roots(1)
      else if (i == 4) then
        z = 1.91d0*z - 0.91d0*positive_roots(2)
      else
        z = 2.d0*z - positive_roots(i-2)
      end if

      do iter = 1, 100
        p1 = PI**(-0.25d0)
        p2 = 0.d0
        do j = 1, n
          p0 = p2
          p2 = p1
          p1 = z*dsqrt(2.d0/dble(j))*p2 &
               - dsqrt(dble(j-1)/dble(j))*p0
        end do
        dp = dsqrt(2.d0*dble(n))*p2
        zold = z
        z = zold - p1/dp
        if (abs(z-zold) <= eps_tol*(abs(z)+1.d0)) exit
      end do
      if (iter > 100) then
        write(*,*) 'gauss_hermite: root iteration did not converge, root ', i
        error stop
      end if
      positive_roots(i) = z
      positive_weights(i) = 2.d0/(dp*dp)

      roots(i) = -z
      roots(n+1-i) = z
      wgh(i) = positive_weights(i)
      wgh(n+1-i) = positive_weights(i)
    end do

    ! store on physical axis: scale nodes by sqrt(Tref); convert weight to
    ! plain ∫ g(v) dv form by multiplying with exp(v^2) and the sqrt(Tref) scale.
    unum = n
    allocate(uspace(unum))
    allocate(weight(unum))
    do i = 1, unum
      uspace(i) = roots(i)*scale
      weight(i) = wgh(i)*dexp(roots(i)**2)*scale
    end do

    ! refresh derived bounds to match the actual node span
    umin = uspace(1)
    umax = uspace(unum)
    du   = (umax-umin)/(unum-1)
    lengthU = umax - umin

    write(*,*) "velocity: Gauss-Hermite, unum =", unum, " Tref =", Tref
    deallocate(roots, wgh, positive_roots, positive_weights)
  end subroutine init_velocity_gauss_hermite
! ------------------------------------------------------------
end module VELOCITY
