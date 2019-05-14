Module postpro
! contains subroutines for post processing
use timeadv

contains

!===================================================================================================

subroutine nusselt_b(T,Nu)
! calculate the Nusselt number of the system
!( Nu = -\int \partial T /\partial y dx ) along bottom plate
implicit none
real*8, dimension(nx-1,ny-1)::T
real*8, dimension(nx-1) :: dTy
real*8 :: Nu

integer :: i,j

dTy(:) = (-3*T(:,1) + 4*T(:,2) - T(:,3))/(2*dy) ! one-side finite difference

call line_int(-dTy,nx-1,Nu)

end subroutine nusselt_b

!===================================================================================================

subroutine nusselt_v1(T,v,Nu)
! (Nu = 1 + sqrt(RaPr)< v T >_V) calculate Nu by volume integral
implicit none
real*8, dimension(nx-1,ny-1)::T,vT
real*8, dimension(nx+1,ny)::v
real*8 :: Nu

integer :: i,j

do i = 1,nx-1
	do j = 1,ny-1
      
		vT(i,j) = (v(i+1,j) + v(i+1,j+1))/2 * T(i,j)

	enddo
enddo

call vol_int(vT,nx-1,ny-1,Nu)

Nu = Nu + 1.0

end subroutine nusselt_v1

!===================================================================================================

subroutine output_Nu(Nu,timet,iid)
implicit none
real*8 :: Nu, timet
integer :: iid

open(iid,file = 'results/Nu_b.dat')
write(iid,*)timet,Nu
end subroutine output_Nu

!===================================================================================================

subroutine line_int(a,n,integ)
! calculate the average of path integral
implicit none
real*8,dimension(n) :: a
real*8 :: integ,len
integer :: i,n

integ = 0.0
len = 0.0

do i = 1,n-1

	integ  = integ + (a(i) + a(i+1))*dx/2 ! trapezoid rule for Riemann integral

	len = len + dx

enddo

integ = integ / len ! average

end subroutine line_int

!===================================================================================================

subroutine vol_int(a,n1,n2,integ)
! calculate the average of volume integral, scalar defined on center of grid cell
implicit none
real*8,dimension(n1,n2) :: a
real*8 :: integ,vol
integer :: i,j,n1,n2


integ = 0.0
vol = 0.0

do i = 1,n1

	do j = 1,n2

		integ  = integ + a(i,j) * dx * dy

		vol = vol + dx * dy

	enddo

enddo

integ = integ/vol ! average

end subroutine vol_int

!===================================================================================================

subroutine convergence(ifcon,u,v,T,u0,v0,T0) ! debug needed
! judge if flow is steady
implicit none
real*8, dimension(nx,ny+1)::u,u0
real*8, dimension(nx+1,ny)::v,v0
real*8, dimension(nx-1,ny-1)::T,T0
real*8 :: err
logical:: ifcon
integer :: i,j

ifcon = .FALSE.

err = 0
if (if_RB ) then
	err = max(err, maxval(abs(u-u0)),maxval(abs(v-v0)),maxval(abs(T-T0)))
else
	err = max(err, maxval(abs(u-u0)),maxval(abs(v-v0)))
endif
write(*,*)err
if (err < 1e-5) then

ifcon = .TRUE.

endif

end subroutine convergence

!===================================================================================================

subroutine recon_p(u,v,p_re)
! reconstruct pressure field from momentum Eqn
implicit none
real*8, dimension(nx,ny+1)::u,hx
real*8, dimension(nx+1,ny)::v,hy
real*8,dimension(nxp,nyp) :: f,p_re

integer :: i,j,n

call advect(u,v,hx,hy)

do i = 2, nxp-1
	do j = 2,nyp-1

		f(i,j) = -((hx(i,j) - hx(i-1,j))/dx + (hy(i,j) - hy(i,j-1))/dy)

	enddo
enddo

if (psolver == 1) then ! SOR
	call solve_sor(p_re,f,n)
elseif (psolver == 2) then ! ADI
	call solve_adi(p_re,f,n)
endif

end subroutine recon_p
end 