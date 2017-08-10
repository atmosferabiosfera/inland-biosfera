#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine mix (xm,tm, x1,t1, x2,t2, x3,t3,npt)
! calorimetrically mixes masses x1,x2,x3 with temperatures
! t1,t2,t3 into combined mass xm with temperature tm
!
! xm,tm may be returned into same location as one of x1,t1,..,
! so hold result temporarily in xtmp,ttmp below
!
! will work if some of x1,x2,x3 have opposite signs, but may 
! give unphysical tm's
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer npt            ! number of points in little vector 

! input-output variables
      real*8 xm(mpt),   & ! resulting mass  
             tm(mpt),   & ! resulting temp
             x1(mpt),   & ! mass 1
             t1(mpt),   & ! temp 1
             x2(mpt),   & ! mass 2
             t2(mpt),   & ! temp 2
             x3(mpt),   & ! mass 3
             t3(mpt)      ! temp 3
!
! local variables
      integer i
      real*8 xtmp, ytmp, ttmp

      do 100 i = 1, npt
        xtmp = x1(i) + x2(i) + x3(i)
        ytmp = sign (max (abs(xtmp), epsilon), xtmp)
        if (abs(xtmp).ge.epsilon) then
          ttmp = (t1(i)*x1(i) + t2(i)*x2(i) + t3(i)*x3(i)) / ytmp
        else
          ttmp = 0.
          xtmp = 0.
        endif
        xm(i) = xtmp
        tm(i) = ttmp
100   continue
      return
end subroutine mix
