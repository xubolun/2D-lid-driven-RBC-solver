Module prms
! parameters used in program
implicit none
integer,parameter :: nx=129 ! grid points number on physical domain 
integer,parameter :: ny=129
integer,parameter :: nxp = nx+1 ! grid points number for p
integer,parameter :: nyp = ny+1  
integer,parameter :: bcflag = 1 ! suggests which type of BCs is imposed 1- wall ; 2- periodic for w/e
integer,parameter :: nstep = 200000 ! total time steps
integer,parameter :: outstep =200! output steps
integer,parameter :: ifre = 0 ! 0 - start from the beginning; others - 'l' for starting from (l)th file
integer,parameter :: chkp = 20000 ! frequency of outputting restart chkpt files
integer,parameter :: torder = 2 ! 1: 1st Euler; 2: 2nd AB
integer,parameter :: psolver = 2 ! 1: SOR; 2: ADI
real*8, parameter :: err = 1e-7 ! error for iteration
real*8, parameter :: dt = 0.0001 ! time step
real*8, parameter :: Ra = 1E4, Pr = 0.71 ! nondimensional parameter
real*8, parameter :: Ut = 1.0, Vt = 0.0 ! lid speed
real*8, parameter :: Ub = 0.0, Vb = 0.0, Ue = 0.0, Ve = 0.0, Uw = 0.0, Vw =0.0 ! other walls, all fixed
real*8, parameter :: Tt = 0.0, Tb = 1.0 ! temperature BC
real*8,parameter :: pi=4.0*atan(1.0)
logical, parameter :: f_abcn = .TRUE. ! FALSE for fully explicit
logical, parameter :: if_RB  = .TRUE. ! FALSE if no energy Eqn considered
real*8 :: Re,RaPr
real*8 :: xmin, xmax, ymin, ymax
real*8 :: dx, dy
real*8, allocatable, dimension(:) :: x,y
contains

subroutine initprms
implicit none

if (if_RB) then
	Re = sqrt(Ra/Pr) ! Re is a responding parameter here
else
	Re = 400 ! can set Re here if no energy Eqn included
endif

RaPr = sqrt(Ra*Pr)

print *, '================   Initialization Starts   ================'
write(*,"('   dt = ',e13.5)"),dt
if (.not.if_RB) then
	write(*,"('   Driven cavity flow')")
	write(*,"('   Re = ',e13.5)"),Re
else
	write(*,"('   RB convection with moving top plate')")
	write(*,"('   Ra = ',e13.5)"),Ra
	write(*,"('   Pr = ',e13.5)"),Pr
	write(*,"('   Speed of top plate: ',e13.5)"),Ut
endif
if (bcflag == 1) then
	write(*,"('   Wall BC for W/E')")
elseif(bcflag == 2) then
	write(*,"('   periodic BC for W/E')")
endif

if (psolver == 1) then
	write(*,"('   SOR method for pressure Poisson Eqn')")
elseif (psolver == 2) then
	write(*,"('   ADI method for pressure Poisson Eqn')")
endif

write(*,"('   converged tolerance for pressure solver: ', e13.5)")err

if (torder == 1) then
	write(*,"('   1st Euler + MAC for momentum Eqn')")
elseif(torder == 2) then

	if(f_abcn) then

		write(*,"('   AB2 + CN + fractional-step for momentum Eqn')")

	else 

		write(*,"('   AB2 explicit + fractional-step for momentum Eqn')")

	endif
endif

print *, '================ Initialization Completed! ================'
end subroutine initprms

end