#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine cascade(kpti, kptj)
! ---------------------------------------------------------------------
! steps intercepted h2o due to drip, precip, and min/max limits
!
! calls steph2o for upper leaves, upper stems and lower veg in
! iurn, adjusting precips at each level
!
! local variables:
! xai(i)   = lai and/or sai for veg component
!            (allows steph2o to work on any veg component)
! rain(i)  = rainfall at appropriate level (modified by steph2o)
! train(i) = temperature of rain (modified by steph2o)
! snow(i)  = snowfall at appropriate level (modified by steph2o)
! tsnow(i) = temperature of snow (modified by steph2o)
      use inland_com1d
      use inland_parameters
      use inland_comveg
      use inland_comatm
      use inland_comcrop, only:xirriga, isimagro

      implicit none

!-----------------------------------------------------------------------
! local variables
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector
      integer npt             ! number of points in little vector

      integer i

      real*8 xai(kpti:kptj),rain(kpti:kptj),train(kpti:kptj),snow(kpti:kptj), &
             tsnow(kpti:kptj),x1(kpti:kptj),x2(kpti:kptj),x3(kpti:kptj), &
             x4(kpti:kptj)

      real*8 twet3, &         ! Function: wet bulb temperature (K)
             twetbulb         ! wet bulb temperature (K)

      npt = kptj - kpti + 1

! adjust rainfall and snowfall rates at above-tree level
!
! set wliqmin, wsnomin -- unlike wliq*max, wsno*max, these are
! part of the lsx numerical method and not from the vegetation
! database, and they are the same for all veg components
!
! the value 0.0010 should be small compared to typical precip rates
! times dtime to allow any intercepted h2o to be initiated, but
! not too small to allow evap rates to reduce wliq*, wsno* to
! that value in a reasonable number of time steps
      wliqmin = 0.0010 * (dtime/3600.) * (wliqumax / 0.2)
      wsnomin = 0.0010 * (dtime/3600.) * (wsnoumax / 2.0)

      do 50 i = kpti, kptj
         rainu(i) = raina(i)
!
! add amount for irrigation  - C. Kucharik 04/11/01
!
        if(isimagro .gt. 0) rainu(i) = rainu(i) + xirriga(i)

! set rain temperature to the wet bulb temperature
         if (ta(i) .gt. tmelt) then
            twetbulb = twet3( ta(i), qa(i), psurf(i) )
         else
            twetbulb = tmelt
         endif
         trainu(i) = max (twetbulb, tmelt)
         x1(i) = 0.0
         x2(i) = max (t12(i), tmelt)
50    continue

      call mix(rainu(kpti),trainu(kpti), rainu(kpti),trainu(kpti), &
               x1(kpti),x2(kpti), vzero, vzero,  npt)

      do 52 i = kpti, kptj
         snowu(i) = snowa(i)
         tsnowu(i) = min (ta(i), tmelt)
         x1(i) = 0.0
         x2(i) = min (t12(i), tmelt)
52    continue

      call mix(snowu(kpti),tsnowu(kpti), snowu(kpti),tsnowu(kpti), &
               x1(kpti),x2(kpti), vzero, vzero, npt)

! set up for upper leaves
      do 100 i = kpti, kptj
         xai(i)   = 2.0 * lai(i,2)
         rain(i)  = rainu(i)
         train(i) = trainu(i)
         snow(i)  = snowu(i)
         tsnow(i) = tsnowu(i)
100   continue

! step upper leaves
      call steph2o(tu(kpti), wliqu(kpti), wsnou(kpti), xai(kpti), &
                   pfluxu(kpti), rain(kpti), train(kpti), &
                   snow(kpti), tsnow(kpti), tdripu, tblowu, &
                   wliqumax, wsnoumax, wliqmin, wsnomin, npt)

! set up for upper stems
! the upper stems get precip as modified by the upper leaves
      do 200 i = kpti, kptj
         xai(i) = 2.0 * sai(i,2)
200   continue

! step upper stems
      call steph2o(ts(kpti),  wliqs(kpti),  wsnos(kpti),  xai(kpti), & 
                   pfluxs(kpti),  rain(kpti), train(kpti), &
                   snow(kpti), tsnow(kpti), tdrips, tblows, &
                   wliqsmax, wsnosmax, wliqmin, wsnomin, npt)

! adjust rainfall and snowfall rates at below-tree level
! allowing for upper-veg interception/drip/belowoff
      do 300 i = kpti, kptj
        x1(i) = fu(i)*rain(i)
        x2(i) = (1.-fu(i))*rainu(i)
        x3(i) = 0.0
        x4(i) = max (t34(i), tmelt)
300   continue

      call mix(rainl(kpti),trainl(kpti), x1(kpti),train(kpti), &
               x2(kpti),trainu(kpti), x3(kpti),x4(kpti), npt)

      do 310 i = kpti, kptj
         x1(i) = fu(i)*snow(i)
         x2(i) = (1.-fu(i))*snowu(i)
         x3(i) = 0.0
         x4(i) = min (t34(i), tmelt)
310   continue
 
      call mix(snowl(kpti),tsnowl(kpti), x1(kpti),tsnow(kpti), &
               x2(kpti),tsnowu(kpti), x3(kpti),x4(kpti),npt)

! set up for lower veg
      do 400 i = kpti, kptj
         xai(i)   = 2.0 * (lai(i,1) + sai(i,1))
         rain(i)  = rainl(i)
         train(i) = trainl(i)
         snow(i)  = snowl(i)
         tsnow(i) = tsnowl(i)
400   continue

! step lower veg
      call steph2o(tl(kpti),  wliql(kpti),  wsnol(kpti),  xai(kpti), &
                   pfluxl(kpti), rain(kpti), train(kpti), snow(kpti), & 
                   tsnow(kpti), tdripl, tblowl, wliqlmax, wsnolmax, &
                   wliqmin, wsnomin,npt)

! adjust rainfall and  snowfall rates at soil level,
! allowing for lower-veg interception/drip/blowoff
      do 500 i= kpti, kptj
         x1(i) = fl(i) * rain(i)
         x2(i) = (1.-fl(i)) * rainl(i)
500   continue

      call mix(raing(kpti),traing(kpti),x1(kpti),train(kpti),x2(kpti), &
               trainl(kpti), vzero,vzero, npt)
      do 510 i = kpti, kptj
         x1(i) = fl(i) * snow(i)
         x2(i) = (1.-fl(i)) * snowl(i)
510   continue

      call mix(snowg(kpti),tsnowg(kpti),x1(kpti),tsnow(kpti),x2(kpti), &
               tsnowl(kpti), vzero,vzero, npt)
      return
end subroutine cascade
