#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine impexp2 (wimp, t, told, iter, kpti, kptj)
! ---------------------------------------------------------------------
!
! sets the implicit vs explicit fraction in turvap calcs for
! seaice or snow skin temperatures, to account for temperatures
! of freezing/melting surfaces being constrained at the melt
! point
!
! unlike impexp, don't have to allow for all h2o 
! vanishing within the timestep
!
! wimp   = implicit fraction (0 to 1) (returned)
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
!-----------------------------------------------------------------------
!
! input variables
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector
      integer iter
      real*8  wimp(lbeg:lend), t(lbeg:lend), told(lbeg:lend)

! local variables
      integer i
      integer npt      ! number of points in little vector

! for first iteration, set wimp to fully implicit, and return
      npt = kptj - kpti + 1
      if (iter .eq. 1) then
         wimp(:)=1.0
         return
      endif

      do 100 i = kpti, kptj
         if ((t(i)-told(i)).gt.epsilon) &
            wimp(i) = (tmelt - told(i)) / (t(i)  - told(i))
         wimp(i) = max (dble(0.0), min (dble(1.0), wimp(i)))
100   continue
      return
end subroutine impexp2
