Module io
! for input & output

USE prms
USE mesh

implicit none

contains

subroutine outputp(p,l)
! output pressure
implicit none
real*8, dimension(nxp,nyp) :: p
real*8, dimension(nx,ny) :: pout
character(200) :: fname_p
integer:: i,j,l

write(fname_p,"('results/p',I7.7,'.q')")l

do j = 1,ny
	do i = 1,nx

		pout(i,j) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4

	enddo
enddo

open(1, file = fname_p,form='UNFORMATTED')
write(1)nx,ny,1
write(1)((pout(i,j),i=1,nx),j=1,ny)
close(1)

end subroutine outputp

subroutine outputT(T,l)
! output temperature
implicit none
real*8, dimension(nx-1,ny-1) :: T
real*8, dimension(nxp,nyp) :: T1
real*8, dimension(nx,ny) :: Tout
character(200) :: fname_T
integer:: i,j,l

T1(2:nx,2:ny) = T(:,:)
T1(1,2:ny) = T(1,:)
T1(nxp,2:ny) = T(nx-1,:)
T1(:,1) = T1(:,2)
T1(:,nyp) = T1(:,nyp-1)

do j = 1,ny
	do i = 1,nx

		Tout(i,j) = (T1(i,j)+T1(i+1,j)+T1(i,j+1)+T1(i+1,j+1))/4

	enddo
enddo

write(fname_T,"('results/T',I7.7,'.q')")l

open(1, file = fname_T,form='UNFORMATTED')
write(1)nx,ny,1
write(1)((Tout(i,j),i=1,nx),j=1,ny)
close(1)
end subroutine outputT

subroutine outputuv(u,v,l)
! output u/v
implicit none
real,dimension(nx,ny) :: uout,vout
real*8, dimension(nx,ny+1) :: u 
real*8, dimension(nx+1,ny) :: v
character(200) :: fname_v
integer:: i,j,l

do j = 1,ny

	uout(:,j) = (u(:,j) + u(:,j+1))/2

enddo

do i = 1,nx

	vout(i,:) = (v(i,:) + v(i+1,:))/2

enddo

write(fname_v,"('results/u',I7.7,'.q')")l

open(1, file = fname_v ,form='UNFORMATTED')
write(1)nx,ny,2
write(1)((uout(i,j),i=1,nx),j=1,ny),((vout(i,j),i=1,nx),j=1,ny)
close(1)
end subroutine outputuv

subroutine checkpoint(u,v,T,timet,l)
! output flow field data for continuing computation
implicit none
real*8, dimension(nx-1,ny-1) :: T
real*8, dimension(nx,ny+1) :: u
real*8, dimension(nx+1,ny) :: v
real*8 :: timet 
character(200) :: fname_re
integer:: l

write(fname_re,"('results/restart/re_',I7.7,'.q')")l
open(1,file = fname_re,form='UNFORMATTED' )
write(1)timet
write(1)u,v,T
close(1)
end subroutine checkpoint

subroutine checkpoint0(u0,v0,T0,timet0,l)
! output flow field data for continuing computation, for 2-step methods
implicit none
real*8, dimension(nx-1,ny-1) :: T0
real*8, dimension(nx,ny+1) :: u0
real*8, dimension(nx+1,ny) :: v0
real*8 :: timet0
character(200) :: fname_re
integer:: l

write(fname_re,"('results/restart/re0_',I7.7,'.q')")l
open(1,file = fname_re,form='UNFORMATTED' )
write(1)timet0
write(1)u0,v0,T0
close(1)
end subroutine checkpoint0

subroutine read_chkpt(u0,v0,T0,u,v,T,timet,l)
! read flow field data from saved files for continuing computation
implicit none
real*8, dimension(nx-1,ny-1) :: T,T0
real*8, dimension(nx,ny+1) :: u,u0
real*8, dimension(nx+1,ny) :: v,v0
real*8 :: timet 
character(200) :: fname_re0, fname_re
integer:: l

write(fname_re0,"('results/restart/re0_',I7.7,'.q')")l
write(fname_re,"('results/restart/re_',I7.7,'.q')")l

open(1,file = fname_re0,form='UNFORMATTED' )
read(1)
read(1)u0,v0,T0
close(1)

open(1,file = fname_re,form='UNFORMATTED' )
read(1)timet
read(1)u,v,T
close(1)
write(*,"('   Restart from ', I7.7, 'th time step...')")l
end subroutine read_chkpt
end