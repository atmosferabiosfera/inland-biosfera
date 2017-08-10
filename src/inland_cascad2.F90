#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine cascad2(kpti, kptj)
! ---------------------------------------------------------------------
!
!        at end of timestep, removes evaporation from intercepted h2o,
!        and does final heat-conserving adjustment for any liquid/snow 
!        below/above melt point. calls steph2o2 for upper leaves, 
!        upper stems and lower veg in turn.
!
!        local variables:
!        (fveg, xai allow steph2o2 to work on any veg component)
!        fveg(i)  = fractional areal coverage of veg component
!        xai(i)   = lai and/or sai for veg component
!
      use inland_com1d
      use inland_parameters
      use inland_comveg
      use inland_comatm
      use inland_comsno

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti  ! index of 1st point of little vector in big lpt vector
      integer kptj  ! index of last point of little vector

! local variables
      integer i
      integer npt             ! number of points little vector

      real*8 fveg(kpti:kptj), &   ! fractional areal coverage of veg component
             xai(kpti:kptj)       ! lai and/or sai for veg component

      npt = kptj - kpti + 1

! set up for upper leaves
      do 100 i= kpti, kptj
         fveg(i) = fu(i)
         xai(i) = 2.0 * lai(i,2)
100   continue

! step upper leaves
      call steph2o2(tu(kpti), wliqu(kpti), wsnou(kpti), fveg(kpti), xai(kpti), &
                    rliqu(kpti), fvapuw(kpti), fvapa(kpti), fsena(kpti), &
                    ta(kpti), chu, npt)

! set up for upper stems
      do 200 i= kpti, kptj
         fveg(i) = fu(i)
         xai(i) = 2.0 * sai(i,2)
200   continue

! step upper stems
      call steph2o2(ts(kpti), wliqs(kpti), wsnos(kpti), fveg(kpti), xai(kpti), &
                    rliqs(kpti), fvaps(kpti), fvapa(kpti), fsena(kpti),        &
                    ta(kpti), chs, npt)

! set up for lower veg
      do 400 i= kpti, kptj
         fveg(i) = (1.-fi(i))*fl(i)
         xai(i) = 2.0 * (lai(i,1) + sai(i,1))
400   continue

! step lower veg
      call steph2o2(tl(kpti), wliql(kpti), wsnol(kpti), fveg(kpti), xai(kpti), &
                    rliql(kpti), fvaplw(kpti), fvapa(kpti), fsena(kpti),       &
                    ta(kpti), chl, npt)
      return
end subroutine cascad2
