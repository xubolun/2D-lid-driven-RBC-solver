Module mesh
! uniform mesh

USE prms

implicit none

contains

subroutine initmesh 
implicit none
integer :: i,j
allocate(x(nx),y(ny))
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0

   dx = (xmax-xmin)/Real(nx-1)   ! calculate dx and dy
   dy = (ymax-ymin)/Real(ny-1)
   
    do i = 1, nx
       x(i) = xmin + dx*Real(i-1)    ! calculate the position of each mesh point
    enddo

    do j = 1, ny
       y(j) = ymin + dy*Real(j-1)
    enddo

end subroutine initmesh

subroutine outmesh
! output grid file
implicit none
integer :: i,j

open(1,file='results/mesh.g',form='UNFORMATTED')
write(1)nx,ny
write(1)((x(i),i=1,nx),j=1,ny),((y(j),i=1,nx),j=1,ny)
close(1)
end subroutine outmesh

end