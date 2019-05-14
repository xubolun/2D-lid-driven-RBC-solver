Module timeadv
!============================================!
! 											 !
! project method for NS time-advancing       !
! containing subroutines for advection term, !
! viscous term, pressure correction.         !
!										     !
!   for energy eqn, AB-CN scheme is used     !
!============================================!
use prms
use mesh
use TDMA
use SOR
use ADI
use BC
use io

implicit none

contains

!===================================================================================================

subroutine ABCN(u,v,u0,v0,T,T0,us,vs,l)
! AB2 + CN scheme, output intermediate velocity us/vs 
! AB2(explicit) for advective + buoyancy(or other body forces) term
! CN(implicit) for viscous term

implicit none
real*8, dimension(nx,ny+1) :: u,u0,hx,hx0,us,fnx,cx,lx,dus,du
real*8, dimension(nx+1,ny) :: v,v0,hy,hy0,vs,fny,cy,ly,dvs,dv,by,by0
real*8, dimension(nx-1,ny-1) :: T,T0
real*8, dimension(nx-2) :: d1u,b1u
real*8, dimension(ny-1) :: d2u,b2u
real*8, dimension(nx-1) :: d1v,b1v
real*8, dimension(ny-2) :: d2v,b2v
real*8 :: a1, a2 

integer :: i,j,l

a1 = dt/(2*Re*dx**2)

a2 = dt/(2*Re*dy**2)

call initBCv(us,vs)

du = 0.0
dv = 0.0
dus = 0.0
dvs = 0.0

if (l>1) then ! if time step>1 then AB2 for advection term

	call advect(u0,v0,hx0,hy0) ! (n-1)th time step

	call advect(u,v,hx,hy) ! nth time step

	call buoyancy(T0,by0) ! (n-1)th time step

	call buoyancy(T,by) ! nth time step

	cx = dt/2*(3*hx-hx0)

	cy = dt/2*((3*hy-hy0) + (3*by-by0)) 

else ! Euler for the first step

	call advect(u,v,hx,hy)

	call buoyancy(T,by)

	cx = dt*hx

	cy = dt*(hy + by)

endif

if (l > 1) then 

	u0 = u ! u_{n} --> u_{n-1}

	v0 = v ! v_{n} --> v_{n-1}

endif

call viscous(u,v,lx,ly) ! nth time step

fnx = cx + 2*dt/2*lx ! rhs vector

fny = cy + 2*dt/2*ly 

! approximate factorization to solve viscous term, use TDMA solver for each direction

! x- direction:
! -a1*dus_{i-1} + (2*a1+1)*dus_{i} - a1*dus_{i+1} = fnx
! -a1*dvs_{i-1} + (2*a1+1)*dvs_{i} - a1*dvs_{i+1} = fny

! construct tridiagonal system

d1u = 2*a1+1.

d1v = 2*a1+1.

! for u
do j = 2,ny

	b1u(1) = fnx(2,j) - (-a1*dus(1,j))
    b1u(nx-2) = fnx(nx-1,j) - (-a1*dus(nx,j))
    b1u(2:nx-3) = fnx(3:nx-2,j)

    call my_thomas(-a1,d1u,b1u,dus(2:nx-1,j),nx-2)

enddo

! for v
do j = 2,ny-1

	b1v(1) = fny(2,j) - (-a1*dvs(1,j))
	b1v(nx-1) = fny(nx,j) - (-a1*dvs(nx+1,j))
	b1v(2:nx-2) = fny(3:nx-1,j)

	call my_thomas(-a1,d1v,b1v,dvs(2:nx,j),nx-1)

enddo


! y- direction:
! -a2*du_{j-1} + (2*a2+1)*du_{j} - a2*du_{j+1} = dus
! -a2*dv_{j-1} + (2*a2+1)*dv_{j} - a2*dv_{j+1} = dvs

! construct tridiagonal system

d2u = 2*a2+1.

d2v = 2*a2+1.

! for u
do i = 2,nx-1

	b2u(1) = dus(i,2) - (-a2*du(i,1))
	b2u(ny-1) = dus(i,ny) - (-a2*du(i,ny+1))
	b2u(2:ny-2) = dus(i,3:ny-1)

	call my_thomas(-a2,d2u,b2u,du(i,2:ny),ny-1)

enddo

! for v
do i = 2,nx

	b2v(1) = dvs(i,2) - (-a2*dv(i,1))
	b2v(ny-2) = dvs(i,ny-1) - (-a2*dv(i,ny))
	b2v(2:ny-3) = dvs(i,3:ny-2)

	call my_thomas(-a2,d2v,b2v,dv(i,2:ny-1),ny-2)

enddo

us = u + du

vs = v + dv

end subroutine ABCN

!===================================================================================================

subroutine AB2exp(u,v,u0,v0,T,T0,us,vs,l)
! AB2 advancing both advective/viscous terms explicitly
implicit none
real*8, dimension(nx,ny+1) :: u,u0,us,hx,lx,hx0,lx0
real*8, dimension(nx+1,ny) :: v,v0,vs,hy,ly,hy0,ly0,by,by0
real*8, dimension(nx-1,ny-1) :: T,T0

integer :: i,j,l


if (l > 1) then ! AB2

	call advect(u0,v0,hx0,hy0) ! (n-1)th time step

	call advect(u,v,hx,hy) ! nth time step

	call viscous(u0,v0,lx0,ly0) ! (n-1)th time step

	call viscous(u,v,lx,ly) ! nth time step

	call buoyancy(T0,by0) ! (n-1)th time step

	call buoyancy(T,by) ! nth time step

	us(2:nx-1,2:ny) = dt/2*(3*(hx(2:nx-1,2:ny)+lx(2:nx-1,2:ny)) - (hx0(2:nx-1,2:ny)+lx0(2:nx-1,2:ny))) + u(2:nx-1,2:ny)

	vs(2:nx,2:ny-1) = dt/2*(3*(hy(2:nx,2:ny-1)+ly(2:nx,2:ny-1)+by(2:nx,2:ny-1)) - (hy0(2:nx,2:ny-1)+ly0(2:nx,2:ny-1)+by0(2:nx,2:ny-1))) + v(2:nx,2:ny-1)

else ! Euler

	call advect(u,v,hx,hy)

	call viscous(u,v,lx,ly)

	call buoyancy(T,by) 

	us(2:nx-1,2:ny) = dt*(hx(2:nx-1,2:ny)+lx(2:nx-1,2:ny)) + u(2:nx-1,2:ny)

	vs(2:nx,2:ny-1) = dt*(hy(2:nx,2:ny-1)+ly(2:nx,2:ny-1)+by(2:nx,2:ny-1)) + v(2:nx,2:ny-1)

endif

if (l > 1) then 

	u0 = u ! u_{n} --> u_{n-1}

	v0 = v ! v_{n} --> v_{n-1}

endif
call initBCv(us,vs)

end subroutine AB2exp

!===================================================================================================

subroutine Euler(u,v,T,us,vs)
! Euler 1st order advancing both advective/viscous terms explicitly
implicit none
real*8, dimension(nx,ny+1) :: u,u0,us,hx,lx
real*8, dimension(nx+1,ny) :: v,v0,vs,hy,ly,by
real*8, dimension(nx-1,ny-1) :: T

integer :: i,j,l


call advect(u,v,hx,hy)

call viscous(u,v,lx,ly)

call buoyancy(T,by)

us(2:nx-1,2:ny) = dt*(hx(2:nx-1,2:ny)+lx(2:nx-1,2:ny)) + u(2:nx-1,2:ny)

vs(2:nx,2:ny-1) = dt*(hy(2:nx,2:ny-1)+ly(2:nx,2:ny-1)+by(2:nx,2:ny-1)) + v(2:nx,2:ny-1)

call initBCv(us,vs)

end subroutine Euler

!===================================================================================================

subroutine advect(u,v,hx,hy)
! discrete advection term hx/hy on staggered mesh
real*8, dimension(nx,ny+1) :: u,hx
real*8, dimension(nx+1,ny) :: v,hy
real*8, dimension(nx-1,ny-1) :: u2 ! u^2
real*8, dimension(nx-1,ny-1) :: v2 ! v^2
real*8, dimension(nx,ny) ::uv ! u*v
integer :: i,j


! hx = -\partial (u^2)/\partial x - \partial (uv)/ \partial y, only compute interior points
! hy = -\partial (uv)/\partial x - \partial (v^2)/ \partial y, only compute interior points

hx = 0.0

hy = 0.0

do i = 1,nx-1

	u2(i,:) = ((u(i,2:ny)+u(i+1,2:ny))/2)**2

enddo

do j = 1, ny-1

	v2(:,j) = ((v(2:nx,j)+v(2:nx,j+1))/2)**2

enddo

do i = 1,nx
	do j = 1,ny

		uv(i,j) = (u(i,j) + u(i,j+1))/2 * (v(i,j) + v(i+1,j))/2

	enddo
enddo

do i = 2,nx-1
	do j = 2,ny

		hx(i,j) = - (u2(i,j-1)-u2(i-1,j-1))/dx - (uv(i,j)-uv(i,j-1))/dy

	enddo
enddo

do i = 2,nx
	do j = 2, ny-1

		hy(i,j) = - (uv(i,j)-uv(i-1,j))/dx - (v2(i-1,j)-v2(i-1,j-1))/dy

	enddo
enddo
end subroutine advect

!===================================================================================================

subroutine viscous(u,v,lx,ly)
! discrete viscous term for nth time step
implicit none
real*8, dimension(nx,ny+1) :: u,lx
real*8, dimension(nx+1,ny) :: v,ly

real*8 :: a1, a2 

integer :: i,j

a1 = 1.0/(Re*dx**2)

a2 = 1.0/(Re*dy**2)

lx = 0.0
ly = 0.0

do i = 2, nx-1
	do j = 2,ny

		lx(i,j) =  a1 * (u(i-1,j) - 2*u(i,j) + u(i+1,j)) + a2 * (u(i,j-1) - 2*u(i,j) + u(i,j+1))

	enddo
enddo

do i = 2, nx
	do j = 2,ny-1

		ly(i,j) =  a1 * (v(i-1,j) - 2*v(i,j) + v(i+1,j)) + a2 * (v(i,j-1) - 2*v(i,j) + v(i,j+1))
		
	enddo
enddo
end subroutine viscous

!===================================================================================================

subroutine pressure(us,vs,p)
! pressure poisson solver
implicit none

real*8,dimension(nxp,nyp) :: f,p
real*8, dimension(nx,ny+1) :: us
real*8, dimension(nx+1,ny) :: vs

integer ::i,j,n

f = 0.0

do i = 2, nxp-1
	do j = 2,nyp-1

		f(i,j) =1/dt*((us(i,j) - us(i-1,j))/dx + (vs(i,j) - vs(i,j-1))/dy)

	enddo
enddo

if (psolver == 1) then ! SOR
	call solve_sor(p,f,n)
elseif (psolver == 2) then ! ADI
	call solve_adi(p,f,n)
endif

end subroutine pressure

!===================================================================================================

subroutine p_correct(us,vs,u,v,p)
! correction on velocity field by gradient of p to get u/v at (n+1) time step
implicit none
real*8, dimension(nx,ny+1) :: u,us
real*8, dimension(nx+1,ny) :: v,vs
real*8, dimension(nxp,nyp) :: p

integer :: i,j

do i = 2,nx-1

	u(i,2:ny) = us(i,2:ny) - dt*(p(i+1,2:ny) - p(i,2:ny))/dx

enddo

do j = 2,ny-1

	v(2:nx,j) = vs(2:nx,j) - dt*(p(2:nx,j+1) - p(2:nx,j))/dy

enddo

end subroutine p_correct

!===================================================================================================

subroutine physical_p(p,preal)
! correction on pressure in flow field (Brown, et al., 2001, JCP)
implicit none
real*8, dimension(nxp,nyp) :: p,preal
integer :: i,j

do i = 2,nxp-1
	do j = 2,nyp-1

		preal(i,j) = p(i,j) - dt/(2*Re)*((p(i+1,j)-2*p(i,j)+p(i-1,j))/dx**2 + (p(i,j+1)-2*p(i,j)+p(i,j-1))/dy**2)

	enddo
enddo

call initBCp(preal)
end subroutine physical_p

!===================================================================================================

subroutine buoyancy(T,by)
! discrete buoyancy force term, in y direction only
implicit none
real*8, dimension(nx+1,ny) :: by
real*8, dimension(nx-1,ny-1) :: T

integer :: i,j

by = 0.0

if (if_RB) then

	do j = 2,ny-1

		by(2:nx,j) = (T(:,j) + T(:,j-1))/2

	enddo

endif
end subroutine buoyancy

!===================================================================================================

subroutine T_adv_ABCN(T1,T,T0,u,v,u0,v0,l)
! temporal advancing energy Eqn with AB/CN scheme
implicit none
real*8, dimension(nx,ny+1) :: u,u0
real*8, dimension(nx+1,ny) :: v,v0
real*8, dimension(nx-1,ny-1) :: T,T0,T1,hT,hT0,lT,c,dtemp,dtemps,fn
real*8, dimension(nx-3) :: d1,b1
real*8, dimension(ny-3) :: d2,b2
real*8 :: a1, a2
integer :: i,j,l

a1 = dt/(2*RaPr*dx**2)

a2 = dt/(2*RaPr*dy**2)

dtemp = 0.0
dtemps = 0.0

if (l>1) then ! if time step>1 then AB2 for advection term

	call advect_T(T0,u0,v0,hT0) ! (n-1)th time step

	call advect_T(T,u,v,hT) ! nth time step

	c = dt/2*(3*hT-hT0)

	T0 = T ! T_{n} --> T_{n-1}

else ! Euler for the first step

	call advect_T(T,u,v,hT)

	c = dt*hT

endif

call diffus_T(T,lT)

fn = c + 2*dt/2*lT ! rhs vector

! approximate factorization to solve diffusive term, use TDMA solver for each direction

! x- direction:
! -a1*dTs_{i-1} + (2*a1+1)*dTs_{i} - a1*dTs_{i+1} = fn

! construct tridiagonal system

d1 = 2*a1+1.

do j = 2,ny-2

	b1(1) = fn(2,j) - (-a1*dtemps(1,j))
    b1(nx-3) = fn(nx-2,j) - (-a1*dtemps(nx-1,j))
    b1(2:nx-4) = fn(3:nx-3,j)

    call my_thomas(-a1,d1,b1,dtemps(2:nx-2,j),nx-3)

enddo

! y- direction:
! -a2*dT_{j-1} + (2*a2+1)*dT_{j} - a2*dT_{j+1} = dTs

! construct tridiagonal system

d2 = 2*a2+1.

! for u
do i = 2,nx-2

	b2(1) = dtemps(i,2) - (-a2*dtemp(i,1))
	b2(ny-3) = dtemps(i,ny-2) - (-a2*dtemp(i,ny-1))
	b2(2:ny-4) = dtemps(i,3:ny-3)

	call my_thomas(-a2,d2,b2,dtemp(i,2:ny-2),ny-3)

enddo

T1 = T + dtemp

call initBCT(T1)

end subroutine T_adv_ABCN

!===================================================================================================

subroutine T_adv_exp(T1,T,T0,u,v,u0,v0,l)
! temporal advancing energy Eqn with AB2 explicit scheme
implicit none
real*8, dimension(nx,ny+1) :: u,u0
real*8, dimension(nx+1,ny) :: v,v0
real*8, dimension(nx-1,ny-1) :: T,T0,hT,hT0,lT,lT0,T1
integer :: i,j,l


if (l>1) then ! if time step>1 then AB2 

	call advect_T(T0,u0,v0,hT0) ! (n-1)th time step

	call advect_T(T,u,v,hT) ! nth time step

	call diffus_T(T0,lT0)

	call diffus_T(T,lT)

	T0 = T

	T1 = dt/2*((3*hT-hT0) + (3*lT-lT0)) + T

	 ! T_{n} --> T_{n-1}

else ! Euler for the first step

	call advect_T(T,u,v,hT)

	call diffus_T(T,lT)

	T1 = dt*(hT + lT) + T

endif

call initBCT(T1)

end subroutine T_adv_exp

!===================================================================================================

subroutine T_adv_Euler(T1,T,u,v)
! temporal advancing energy Eqn with 1st order Euler scheme
implicit none
real*8, dimension(nx,ny+1) :: u
real*8, dimension(nx+1,ny) :: v
real*8, dimension(nx-1,ny-1) :: T,hT,lT,T1

call advect_T(T,u,v,hT)

call diffus_T(T,lT)

T1 = dt*(hT + lT) + T

call initBCT(T1)

end subroutine T_adv_Euler

!===================================================================================================

subroutine advect_T(T,u,v,hT)
! discrete advective term in energy Eqn
implicit none
real*8, dimension(nx,ny+1) :: u
real*8, dimension(nx+1,ny) :: v
real*8, dimension(nx-1,ny-1) :: T,hT
real*8, dimension(nx,ny-1) :: uT
real*8, dimension(nx-1,ny) :: vT

integer :: i,j

hT = 0.0

do i = 2,nx-1

        uT(i,2:ny-2) = u(i,3:ny-1)*(T(i,2:ny-2) + T(i-1,2:ny-2))/2
enddo

do j = 2,ny-1

        vT(2:nx-2,j) = v(3:nx-1,j)*(T(2:nx-2,j) + T(2:nx-2,j-1))/2
	
enddo


do j = 2,ny-2
	do i = 2,nx-2

		hT(i,j) = -((uT(i+1,j) - uT(i,j))/dx + (vT(i,j+1) - vT(i,j))/dy)

	enddo
enddo

end subroutine advect_T
!===================================================================================================

subroutine diffus_T(T,lT)
! discrete diffusive term in energy Eqn
implicit none
real*8, dimension(nx-1,ny-1) :: T,lT

integer :: i,j

lT = 0.0

do i = 2,nx-2
	do j = 2,ny-2

		lT(i,j) = 1/RaPr*((T(i-1,j)-2*T(i,j)+T(i+1,j))/dx**2 + (T(i,j-1)-2*T(i,j)+T(i,j+1))/dy**2)

	enddo
enddo
end subroutine diffus_T
!===================================================================================================
end