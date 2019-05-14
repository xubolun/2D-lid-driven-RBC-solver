Module SOR
! solve by SOR method (based on GS)
USE prms
USE mesh
USE BC

implicit none

contains

subroutine solve_sor(phi,del,n_opt)
implicit none
real*8 :: beta,error,rho,omega_opt
real*8, dimension(nxp,nyp) :: phi0,phi,del
integer :: i,j,n_opt
    
beta = dx/dy

 
rho = 0.5* (cos(pi/nxp)+cos(pi/nyp)) ! spectral radius

omega_opt = 2.0/(1.0+sqrt(1.0 - rho**2)) ! theoretical optimal omega

n_opt = 0

! run with theoretical optimal omega

error = 1.0

phi0 = phi

do while (error > err)

  n_opt = n_opt + 1

  error = 0.0

  do j = 2,nyp-1
    do i = 2,nxp-1

      phi(i,j) =(1.0-omega_opt)*phi0(i,j) + omega_opt*(phi0(i+1,j) + phi(i-1,j) + beta**2 *(phi0(i,j+1)+phi(i,j-1)) -dx**2*del(i,j))/(2.0*(1.0+beta**2))

      error = max(error,abs(phi(i,j)-phi0(i,j)))

      phi0(i,j) = phi(i,j)
    enddo
  enddo

  call initBCp(phi)

if(n_opt>1e5)then
  write(*,*)'too many iterations, diverged!'
  exit
endif
enddo

end subroutine solve_sor

end