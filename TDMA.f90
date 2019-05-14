Module TDMA
! Thomas algorithm for solving tri-diagonal system
implicit none

contains

subroutine my_thomas(a,b,c,x,n) ! for certain case only
! input :
! a: coefficients in upper & lower diagonal, same values
! b: coefficients in main diagonal
! c: rhs vector
! n: number of equations
! output:
! x(n)   - solutions
implicit none 
integer :: n
real*8 :: a
real*8, dimension(n) :: c,b,x
integer :: i

!step 1: forward elimination
do i=2,n
  b(i) = b(i) - (a * a) / b(i-1)
  c(i) = c(i) - c(i-1) * a / b(i-1)
enddo

!step 2: backward substitution
x(n) = c(n)/b(n)
do i=n-1,1,-1
   x(i) = (c(i)- a * x(i+1))/b(i)
enddo
end subroutine my_thomas
end