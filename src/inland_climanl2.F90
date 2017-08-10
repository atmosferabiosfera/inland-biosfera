#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
subroutine climanl2
! ---------------------------------------------------------------------
!
! this subroutine updates the growing degree days, coldest temp, and
! warmest temp if monthly anomalies or daily values are used
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comveg
      use inland_combcs

      implicit none
! ---------------------------------------------------------------------
!
! input variables
!
      ! integer isimveg (removed)
!
! local variables
!
      integer i
      real*8 zweigc, zweigw, rworkc, rworkw
!
!     The filtering of the growing degree days and the climatic limits
!     for existence of pft is done over 5 years instead of 30 in off
!     -line inland.
!
      zweigc = exp(-1./5.)
      zweigw = exp(-1./5.)
!
      rworkc = 1. - zweigc
      rworkw = 1. - zweigw
!
! update critical climatic parameters with running average
!
      do 100 i = lbeg, lend
!     
         tc(i) = zweigc * tc(i) + rworkc * tcthis(i)
         tw(i) = zweigw * tw(i) + rworkw * twthis(i)
!
         tcmin(i) = tc(i) + deltat(i)
!
         gdd0(i) = zweigc * gdd0(i) + rworkc * gdd0this(i)
!
         gdd5(i) = zweigc * gdd5(i) + rworkc * gdd5this(i)
!
! Initialize this year's value of gdd0, gdd5, tc and tw to 0
! (climanl2 called 1st time step of the year, different from offline inland)
!
         tcthis(i) = 100.
         twthis(i) = - 100.
         gdd0this(i) = 0.
         gdd5this(i) = 0.
!
 100  continue
!
      call existence
!
      return
end subroutine climanl2
