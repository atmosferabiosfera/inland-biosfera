#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine tridiaib(ns, nd, ne, a, b, c, y, x, alpha, gamma)
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
!-----------------------------------------------------------------------
!     purpose:
!     to compute the solution of many tridiagonal linear systems.
!
!      arguments:
!      ns ..... the number of systems to be solved.
!
!      nd ..... first dimension of arrays (ge ns) - unused.
!
!      ne ..... the number of unknowns in each system.
!               this must be > 2. second dimension of arrays.
!
!      a ...... the subdiagonals of the matrices are stored
!               in locations a(j,2) through a(j,ne).
!
!      b ...... the main diagonals of the matrices are stored
!               in locations b(j,1) through b(j,ne).
!
!      c ...... the super-diagonals of the matrices are stored in
!               locations c(j,1) through c(j,ne-1).
!
!      y ...... the right hand side of the equations is stored in
!               y(j,1) through y(j,ne).
!
!      x ...... the solutions of the systems are returned in
!               locations x(j,1) through x(j,ne).
!
!      alpha .. work array dimensioned alpha(nd,ne)
!
!      gamma .. work array dimensioned gamma(nd,ne)
!
!       history:  based on a streamlined version of the old ncar
!                 ulib subr trdi used in the phoenix climate
!                 model of schneider and thompson (j.g.r., 1981).
!                 revised by starley thompson to solve multiple
!                 systems and vectorize well on the cray-1.
!                 later revised to include a parameter statement
!                 to define loop limits and thus enable cray short
!                 vector loops.
!
!       algorithm:  lu decomposition followed by solution.
!                   note: this subr executes satisfactorily
!                   if the input matrix is diagonally dominant
!                   and non-singular.  the diagonal elements are
!                   used to pivot, and no tests are made to determine
!                   singularity. if a singular or numerically singular
!                   matrix is used as input a divide by zero or
!                   floating point overflow will result.
!
!       last revision date:      4 february 1988
!
! input- output variables
      integer ns, nd, ne

!     real*8 a(nd,ne), b(nd,ne), c(nd,ne), y(nd,ne),&
!          x(nd,ne), alpha(nd,ne), gamma(nd,ne)
      real*8 a(mpt,ne), b(mpt,ne), c(mpt,ne), y(mpt,ne), x(mpt,ne), &
             alpha(mpt,ne), gamma(mpt,ne)

! local variables
      integer nm1, j, i, ib
      nm1 = ne-1

! obtain the lu decompositions
      do 10 j=1,ns
         alpha(j,1) = 1./b(j,1)
         gamma(j,1) = c(j,1)*alpha(j,1)
10    continue
      do 11 i=2,nm1
         do 12 j=1,ns
            alpha(j,i) = 1./(b(j,i)-a(j,i)*gamma(j,i-1))
            gamma(j,i) = c(j,i)*alpha(j,i)
12       continue
11    continue

! solve
      do 20 j=1,ns
         x(j,1) = y(j,1)*alpha(j,1)
20    continue
      do 21 i=2,nm1
         do 22 j=1,ns
            x(j,i) = (y(j,i)-a(j,i)*x(j,i-1))*alpha(j,i)
22       continue
21    continue
      do 23 j=1,ns
         x(j,ne) = (y(j,ne)-a(j,ne)*x(j,nm1)) / (b(j,ne)-a(j,ne)*gamma(j,nm1))
23    continue
      do 24 i=1,nm1
         ib = ne-i
         do 25 j=1,ns
            x(j,ib) = x(j,ib)-gamma(j,ib)*x(j,ib+1)
25       continue
24    continue
      return
end subroutine tridiaib
