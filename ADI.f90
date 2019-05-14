Module ADI
! solve by ADI method
USE prms
USE mesh
USE BC
USE TDMA

implicit none

contains


subroutine solve_adi(phif,del,n) ! based on SLOR
implicit none
real*8 :: beta,error,omega,Ax,Ay
real*8, dimension(nxp,nyp) :: phi0,phi,del,phif
real*8, dimension(nxp-2) :: Bx,Cx
real*8, dimension(nyp-2) :: By,Cy
real*8, dimension((nxp-2)*(nyp-2))  :: R
integer :: i,j,n,k

call initBCp(phi)

call initBCp(phif)

beta = dx/dy

error = 1.0

phi0 = phi

omega = 1.13

n = 0

! construct tri-diagonal matrix

Ax = - omega / (2.0 * (1.0 + beta**2)) ! upper & lower diagonals

Ay = - omega * beta**2 / (2.0 * (1.0 + beta**2))

do while (error>=err)
   k=1
   error = 0.
      do j = 2,nyp-1 ! sweep by row, from j = 2 to j = nyp-1

         Bx= 1. ! main diagonal

         do i = 2,nxp-1 ! line relaxation
         Cx(i-1) = (1.0-omega) * phi(i, j) + omega * 1.0/(2.0 * (1.0+beta**2))*(beta**2 * (phi0(i, j+1) + phi(i, j-1))- dx**2*del(i,j)) 
         enddo

         Cx(1) = Cx(1) - Ax * phi(1,j)

         Cx(nxp-2) = Cx(nxp-2) - Ax * phi(nxp, j)


         call my_thomas(Ax,Bx,Cx,phi(2:nxp-1,j),nxp-2) ! solving tri-diagonal system

      enddo ! get intermediate phi
      
      do i = 2,nxp-1 ! sweep by column, from i = 2 to i = nxp-1

         By = 1. ! main diagonal

         do j = 2,nyp-1 ! line relaxation
         Cy(j-1) = (1.0-omega) * phif(i, j) + omega * 1.0/(2.0 * (1.0+beta**2))*( (phi(i+1, j) + phif(i-1, j))- dx**2*del(i,j)) 
         enddo

         Cy(1) = Cy(1) - Ay * phif(i,1)

         
         Cy(nyp-2) = Cy(nyp-2) - Ay * phif(i, nyp)


         call my_thomas(Ay,By,Cy,phif(i,2:nyp-1),nyp-2) ! solving tri-diagonal system

         do j = 2,nyp-1
         R(k) = phif(i,j) - phi0(i,j) ! residual
         k = k+1
         phi0(i,j) = phif(i,j)
         enddo
     enddo

   
      call initBCp(phif)

      error = max(error,MAXVAL(abs(R)))
   if(n>1e5)then
   write(*,*)'too many iterations, diverged!'
   exit
   endif
   n=n+1
 enddo

end subroutine solve_adi
end