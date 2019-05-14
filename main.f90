program main
!=============================================================!
!                                                             !
!                                                             !
!                      ME820 project2                         !
!                                                             !
!      2D Lid-driven cavity with temperature difference       !
!                                                             !
!             2nd order AB for convection term                !
!             2nd order CN for viscous term                   !
!                      2nd FD method                          !
!                                                             !
!                    2 types of x- BC:                        !
!                          									  !
!    	     1: wall        2: periodic (to be added)         !
!                                                             !
!															  !
! 					by Bolun Xu @ KSU                         !
!															  !
!				        Apr. 2019                             !
!															  !
!=============================================================!
!                                                             !
!                     updating history                        !
!                                                             !
! V1.0 including NS solver for thermal convetion       05/09  !
!=============================================================!
use prms
use mesh
use BC
use io
use trace
use timeadv
use postpro

implicit none

!   staggered mesh used
!
!    ------ v ------
!     |           |
!     |           |
!     u    p/T    u
!     |           |
!     |           |
!    ------ v ------

real*8, dimension(nxp,nyp) :: p,preal ! pressure on grid center
real*8, dimension(nx-1,ny-1) :: T,T0,T1
real*8, dimension(nx,ny+1) :: u0,u,us! center of y grid
real*8, dimension(nx+1,ny) :: v0,v,vs ! center of x grid
real*8 :: timet ! total time in problem
real*8 :: Nu ! Nusselt number
real*8, allocatable, dimension(:) :: x0,y0 ! coordinates of tracing pts
integer :: tstart,tend,clock_rate, clock_max
integer :: l,n 
logical :: ifcon ! if steady state reached
! ======initial flow field==========
ifcon = .FALSE.
call system_clock ( tstart, clock_rate, clock_max )
call initprms
call initmesh 
call outmesh
call init_tracing(x0,y0,n)

if (ifre == 0) then

	u = 0.0
	v = 0.0
	u0 = u
	v0 = v 
	timet = 0.
	call initT(T)
!	T = 0.0
	T0 = T

else

	call read_chkpt(u0,v0,T0,u,v,T,timet,ifre)

endif
p = 0.0

l = 0

write(*,"('processing:  ')", advance = "no")

 call trace_pts(u,v,p,T,x0,y0,n,timet)

do l = 1,nstep

	call initBCv(u,v) ! enforce BC in every iteration
	call initBCT(T)


	if (mod(l,chkp) == 0 .and. torder == 2) then ! save flow field at (n-1)th step for restart, 2-step method only

		call checkpoint0(u0,v0,T0,timet,l)

	endif

if (if_RB) then ! temporal advancing of temperature

	if (f_abcn) then

		call T_adv_ABCN(T1,T,T0,u,v,u0,v0,l)

	else
		if (torder == 1) then

			call T_adv_Euler(T1,T,u,v)

		elseif (torder == 2) then

			call T_adv_exp(T1,T,T0,u,v,u0,v0,l)

		endif

	endif

    T = T1

endif

if (torder == 1) then ! temporal advancing of velocity

	call Euler(u,v,T,us,vs) ! 1st order Euler

elseif (torder == 2) then

    if (f_abcn) then

		call ABCN(u,v,u0,v0,T,T0,us,vs,l) ! AB2 for nonlinear term, CN for linear term, output intermediate velocity

	else

   		call AB2exp(u,v,u0,v0,T,T0,us,vs,l) ! AB2 for both nonlinear/linear terms, output intermediate velocity

   	endif

endif
	
    call pressure(us,vs,p) ! solving pressure 

	call p_correct(us,vs,u,v,p) ! update u/v at (n+1) time step

	call physical_p(p,preal) ! real pressure

	if (mod(l,outstep) == 0) then ! output 

		call outputp(preal,l)
		call outputT(T,l)
		call outputuv(u,v,l)

	endif

	timet = timet + dt

	if (mod(l,chkp) == 0) then ! save flow field at nth step for restart

		call checkpoint(u,v,T,timet,l)

	endif

	 call trace_pts(u,v,p,T,x0,y0,n,timet)

!   ================ Nusselt number =================
	if (if_RB) then
 
		call nusselt_b(T,Nu) ! line-integral over bottom plate 

		call output_Nu(Nu,timet,10)

	endif

!   =================================================

	if (floor(real(25*l/nstep))>floor(real(25*(l-1)/nstep))) write(*,"('>')", advance = "no")

!	call convergence(ifcon,u,v,T,u0,v0,T0)

	if (ifcon) then
		write(*,"('   Steady state reached! ')")
		exit
	endif

enddo

close(10)

! ==================================
print *, ''
print *, '================  Computation Completed!   ================'
call system_clock ( tend, clock_rate, clock_max )
write(*,"('   Elapsed total time in problem: ', e13.5)")timet
write(*,"('   Elapsed total Wall time: ', e13.5,' secs')")real( tend - tstart )/real( clock_rate )
end program main
