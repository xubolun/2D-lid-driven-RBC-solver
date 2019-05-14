Module BC
! initialize BC
USE prms
USE mesh
use ifport

implicit none


contains

subroutine initBCv(u,v)
! BC for velocity
implicit none
real*8, dimension(nx,ny+1) :: u ! center of y grid
real*8, dimension(nx+1,ny) :: v ! center of x grid

if (bcflag == 1) then ! wall BC 

	u(1,:) = Uw ! west BC

	u(nx,:) = Ue ! east BC

	v(1,:) = 2*Vw - v(2,:) ! west BC

	v(nx+1,:) = 2*Ve - v(nx,:) ! east BC

elseif(bcflag == 2) then ! periodic in x-

	u(nx,:) = u(1,:) 

	v(nx+1,:) = v(1,:) 

	
endif

u(:,ny+1) = 2*Ut - u(:,ny)
u(:,1) = 2*Ub - u(:,2)
v(:,ny) = Vt
v(:,1) = Vb

end subroutine initBCv

subroutine initBCp(p)
! BC for pressure
implicit none
real*8, dimension(nxp,nyp) :: p 

! Neumann BC for p: \partial p/ \partial n = 0
if  (bcflag == 1) then
	p(1,:) = p(2,:)
	p(nxp,:) = p(nxp-1,:)
elseif (bcflag == 2) then
	p(nxp,:) = p(1,:)
endif

p(:,1) = p(:,2)
p(:,nyp) = p(:,nyp-1)

end subroutine initBCp

subroutine initBCT(T)
! BC for temperature
real*8, dimension(nx-1,ny-1) :: T

 if (bcflag == 1) then ! adiabatic BC for W/E: \partial T/ \partial n = 0
	T(1,:) = T(2,:)
	T(nx-1,:) = T(nx-2,:)
elseif (bcflag == 2) then ! periodic in x-
	T(nx-1,:) = T(1,:)
endif

! wall BC for top/bottom
T(:,ny-1) = Tt
T(:,1) = Tb

! if (bcflag == 1) then ! adiabatic BC for W/E: \partial T/ \partial n = 0
! 	T(:,1) = T(:,2)
! 	T(:,ny) = T(:,ny-1)
! elseif (bcflag == 2) then ! periodic in x-
! 	T(:,ny) = T(:,1)
! endif
! ! wall BC for l/r, for vertical cavity
! T(1,:) = Tt
! T(nx,:) = Tb

end subroutine initBCT

!===================================================================================================

subroutine initT(T)
! Initialize temperature field at beginning, linear distribution of T + sinusoidal perturbation
implicit none
real*8, dimension(nx-1,ny-1) :: T ! at grid point

integer :: i,j

do j = 1,ny-1

	T(:,j) = (Tt-Tb)/(ny-2)*j+ Tb  - (Tt-Tb)/(ny-2)

enddo

do i = 1,nx-1

	T(i,:) = T(i,:) + random(0)*0.005*sin(2*pi*(x(i)+x(i+1))/2)

enddo
end subroutine initT

!===================================================================================================
end