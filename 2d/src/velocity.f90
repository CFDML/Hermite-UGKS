! ------------------------------------------------------------
! module VELOCITY
! 2D velocity-space node/weight generation for the UGKS solver.
! Two quadratures, selected by quad_type from module CONFIG:
!   'newton_cotes'  : uniform tensor-product nodes, 6th-order closed NC weights
!   'gauss_hermite' : tensor product of 1D Gauss-Hermite rules
! Both produce (uspace, vspace, weight) on the physical velocity axes,
! with weights as plain ∫∫ g(u,v) du dv quadrature weights.
! ------------------------------------------------------------
module VELOCITY
use ALLDATA
use CONFIG
implicit none

contains
! ------------------------------------------------------------
  subroutine init_velocity()
    select case(trim(quad_type))
    case('newton_cotes')
      call init_velocity_newton_2d()
    case('gauss_hermite')
      call init_velocity_gauss_hermite_2d()
    case default
      write(*,*) "unknown quad_type '",trim(quad_type),"'; using newton_cotes"
      call init_velocity_newton_2d()
    end select
  end subroutine init_velocity
! ------------------------------------------------------------
  subroutine init_velocity_newton_2d()
    integer :: i, j
    unum = nv
    vnum = nv
    allocate(uspace(unum,vnum))
    allocate(vspace(unum,vnum))
    allocate(weight(unum,vnum))

    du = (umax-umin)/(unum-1)
    dv = (vmax_-vmin)/(vnum-1)
    do j = 1, vnum
      do i = 1, unum
        uspace(i,j) = umin + (i-1)*du
        vspace(i,j) = vmin + (j-1)*dv
        weight(i,j) = (newton_coeff(i,unum)*du)*(newton_coeff(j,vnum)*dv)
      end do
    end do
    write(*,*) "velocity: 2D Newton-Cotes, unum=vnum =", unum, " du =", du
  end subroutine init_velocity_newton_2d
! ------------------------------------------------------------
  pure function newton_coeff(idx,num)
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
  subroutine init_velocity_gauss_hermite_2d()
    ! Tensor product of the verified 1D physicist-Hermite rule.  Roots are
    ! generated from the outside inward so every Newton solve converges to a
    ! distinct root.  The normalized recurrence avoids overflow.
    integer :: n, i, j, iter, m
    real(8) :: z, zold, p0, p1, p2, dp, scale
    real(8), allocatable :: roots(:), wgh(:)
    real(8), allocatable :: positive_roots(:), positive_weights(:)
    real(8), parameter :: eps_tol = 1.d-14

    n = nv
    if (n < 2) then
      write(*,*) "gauss_hermite: nv must be >= 2, got ", n
      stop
    end if
    m = (n+1)/2
    allocate(roots(n), wgh(n), positive_roots(m), positive_weights(m))
    scale = dsqrt(Tref)

    z = 0.d0
    do i = 1, m
      if (i == 1) then
        z = dsqrt(2.d0*n+1.d0)-1.85575d0*(2.d0*n+1.d0)**(-1.d0/6.d0)
      else if (i == 2) then
        z = z-1.14d0*dble(n)**0.426d0/z
      else if (i == 3) then
        z = 1.86d0*z-0.86d0*positive_roots(1)
      else if (i == 4) then
        z = 1.91d0*z-0.91d0*positive_roots(2)
      else
        z = 2.d0*z-positive_roots(i-2)
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
        z = zold-p1/dp
        if (abs(z-zold) <= eps_tol*(abs(z)+1.d0)) exit
      end do
      if (iter > 100) error stop 'gauss_hermite: root iteration did not converge'
      positive_roots(i) = z
      positive_weights(i) = 2.d0/(dp*dp)
      roots(i) = -z
      roots(n+1-i) = z
      wgh(i) = positive_weights(i)
      wgh(n+1-i) = positive_weights(i)
    end do

    unum = n
    vnum = n
    allocate(uspace(unum,vnum))
    allocate(vspace(unum,vnum))
    allocate(weight(unum,vnum))
    do j = 1, vnum
      do i = 1, unum
        uspace(i,j) = roots(i)*scale
        vspace(i,j) = roots(j)*scale
        weight(i,j) = wgh(i)*dexp(roots(i)**2)*scale * wgh(j)*dexp(roots(j)**2)*scale
      end do
    end do
    umin = uspace(1,1); umax = uspace(unum,1)
    vmin = vspace(1,1); vmax_ = vspace(1,vnum)
    du = (umax-umin)/(unum-1); dv = (vmax_-vmin)/(vnum-1)
    lengthU = umax-umin; lengthV = vmax_-vmin
    write(*,*) "velocity: 2D Gauss-Hermite, unum=vnum =", unum, " Tref =", Tref
    deallocate(roots, wgh, positive_roots, positive_weights)
  end subroutine init_velocity_gauss_hermite_2d
! ------------------------------------------------------------
end module VELOCITY
