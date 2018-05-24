#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine linsolve(arr, rhs, vec, mplate, nd, npt)
! ---------------------------------------------------------------------
!
! solves multiple linear systems of equations, vectorizing
! over the number of systems. basi! gaussian elimination is
! used, with no pivoting (relies on all diagonal elements
! being and staying significantly non-zero)
!
! a template array mplate is used to detect when an operation
! is not necessary (element already zero or would add zeros),
! assuming that every system has the same pattern of zero
! elements
!
! this template is first copied to mplatex since it
! must be updated during the procedure in case an original-zero
! pattern location becomes non-zero
!
! the first subscript in arr, rhs, ve! is over the multiple
! systems, and the others are the usual row, column subscripts
!
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
! -----------------------------------------------------------------------
! input variables
      integer :: npt                 ! number of points in little vector
      integer :: nd                  ! number of equations (supplied)
      integer ::  mplate(nd,nd)      ! pattern of zero elements of arr (supplied)
      real*8  arr(npt,nd,nd), &     ! equation coefficients (supplied, overwritten)
              rhs(npt,nd),    &     ! equation right-hand sides (supplied, overwritten)
              vec(npt,nd)           ! solution (returned)

! local variables
      integer :: ndx, &               ! Max number of equations
                 j, i, id, m          ! loop indices
      parameter (ndx=9)
      real*8 :: f(npt)
      integer :: mplatex(ndx,ndx)
      if (nd.gt.ndx) then
         write(*,900) nd, ndx
900      format(/' *** fatal error *** number of linsolve eqns',i4, &
                ' exceeds limit',i4)
         stop
      endif

! copy the zero template so it can be changed below
      do 6 j=1,nd
         do 5 i=1,nd
           mplatex(i,j) = mplate(i,j)
5       continue
6    continue

! zero all array elements below the diagonal, proceeding from
! the first row to the last. note that mplatex is set non-zero
! for changed (i,j) locations, in loop 20
     do 10 id=1, nd-1
        do 12 i=id+1,nd
           if (mplatex(i,id).ne.0) then
              do 14 m=1, npt
                 f(m) = arr(m,i,id) / arr(m,id,id)
14            continue
              do 20 j=id,nd
                 if (mplatex(id,j).ne.0) then
                    do 22 m=1, npt
                       arr(m,i,j) = arr(m,i,j) - f(m)*arr(m,id,j)
22                  continue
                    mplatex(i,j) = 1
                 endif
20            continue
              do 30 m=1, npt
                 rhs(m,i) = rhs(m,i) - f(m)*rhs(m,id)
30            continue
           endif
12      continue
10   continue

! all array elements below the diagonal are zero, so can
! immediately solve the equations in reverse order
     do 50 id=nd,1,-1
        f(:)=0.
        if (id.lt.nd) then
           do 52 j=id+1,nd
              if (mplatex(id,j).ne.0) then
                 do 54 m=1,npt
                    f(m) = f(m) + arr(m,id,j)*vec(m,j)
54               continue
              endif
52         continue
        endif
        do 56 m=1,npt
           vec(m,id) = (rhs(m,id) - f(m)) / arr(m,id,id)
56      continue
50  continue
    return
end subroutine linsolve
