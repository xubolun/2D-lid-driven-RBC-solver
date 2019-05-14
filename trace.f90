Module trace
use prms
use mesh
implicit none

contains
!===================================================================================================
subroutine init_tracing(x0,y0,n)
real*8, allocatable, dimension(:) :: x0,y0 ! coordinates of tracing pts
integer :: n ! number of pts
integer :: i
integer :: ios
open(1,file = 'pts.in',iostat = ios)
n = 0
! count the number of pts
do
    read(1, * , iostat=ios) 

    if (ios /= 0) exit

    n = n + 1
enddo

write(*,"('   ',I, ' tracing pts in total')"),n

allocate(x0(n),y0(n))

rewind(1)

do i = 1,n
	read(1,*)x0(i),y0(i) 
enddo

close(1)

end subroutine init_tracing
!===================================================================================================

subroutine trace_pts(u,v,p,T,x0,y0,n,timet)
! trace the temporal series of u/v/p/T at given points
implicit none
real*8,dimension(nx,ny+1)  :: u
real*8,dimension(nx+1,ny)  :: v
real*8,dimension(nx-1,ny-1) :: T
real*8,dimension(nxp,nyp) :: p,T1
real*8,dimension(nx,ny)  :: up,vp,pp,Tp
real*8,dimension(n) :: x0,y0
real*8,dimension(4) :: val0
real*8 :: timet
integer :: i,j,l,n 
logical :: exist
character(80), allocatable, dimension(:) :: filename

allocate(filename(n))

do j = 1,ny

	up(:,j) = (u(:,j) + u(:,j+1))/2

enddo

do i = 1,nx

	vp(i,:) = (v(i,:) + v(i+1,:))/2

enddo

do j = 1,ny
	do i = 1,nx

		pp(i,j) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4

	enddo
enddo

T1(2:nx,2:ny) = T(:,:)
T1(1,2:ny) = T(1,:)
T1(nxp,2:ny) = T(nx-1,:)
T1(:,1) = T1(:,2)
T1(:,nyp) = T1(:,nyp-1)

do j = 1,ny
    do i = 1,nx

        Tp(i,j) = (T1(i,j)+T1(i+1,j)+T1(i,j+1)+T1(i+1,j+1))/4

    enddo
enddo

do l = 1,n
	call intp_pts(up,val0(1),x0(l),y0(l))
	call intp_pts(vp,val0(2),x0(l),y0(l))
	call intp_pts(pp,val0(3),x0(l),y0(l))
	call intp_pts(Tp,val0(4),x0(l),y0(l))
    write(filename(l),"('results/pts_', I2.2, '.dat')")l
    inquire(file = filename(l), exist = exist)

    if(exist) then
    	open(1,file = filename(l),status="old", position="append", action="write")
    else
    	open(1,file = filename(l),status="new", action="write")
    endif

    write(1,"(5e13.5)")timet,(val0(i),i=1,4)

    close(1)

enddo

end subroutine trace_pts

!===================================================================================================

subroutine intp_pts(val0,val,x0,y0)
! interpolate flow field on certain pts
implicit none
real*8, dimension(nx,ny) :: val0
real*8 :: x0,y0,val

integer :: i,j

do i = 1,nx-1
	do j = 1,ny-1
		if(x0 >= x(i) .and. x0 <= x(i+1) .and. y0 >= y(j) .and. y0 <= y(j+1) ) then
		!--------------------------------------------------------------------------------------------------
	    ! bilinear interpolation
        ! val = 1/((x(i+1)-x(i))(y(i+1)-y(i))) [x(i+1)-x0 x0-x(i)] [val0(i,j)    val0(i,j+1)  ] [y(i+1)-y0]
        !                                                          [val0(i+1,j)  val0(i+1,j+1)]	[y0 - y(i)]
        !--------------------------------------------------------------------------------------------------

        val = 1/(dx*dy) *( ( (x(i+1)-x0)*val0(i,j) + (x0-x(i))*val0(i+1,j) ) * (y(j+1)-y0)  &
             + ( (x(i+1)-x0)*val0(i,j+1) + (x0-x(i))*val0(i+1,j+1) ) * (y0 - y(j)) )               
		endif
	enddo
enddo
end subroutine intp_pts
!===================================================================================================
end