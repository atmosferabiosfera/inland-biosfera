#include "inland_config.h"
! --------------------------------------------------------------------------
subroutine steph2o(tveg,wliq,wsno,xai,pflux,rain,train,snow,tsnow,tdrip,tblow, &
                   wliqmax,wsnomax,wliqmin,wsnomin,npt)
! --------------------------------------------------------------------------
!
! steps intercepted h2o for one canopy component (upper leaves, 
! upper stems, or lower veg) through one lsx time step, adjusting
! for h2o sensible heat and phase changes. also modifies precip
! due to interception and drip,blowoff
      use inland_parameters

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer npt            ! number of points in little vector 

!  Arguments (all arguments are supplied (unchanged) unless otherwise noted
      real*8 tdrip,   &  ! e-folding time of liquid drip  tdrip[u,s,l]
             tblow,   &  ! e-folding time of snow blowoff tblow[u,s,l]
             wliqmax, &  ! max amount of intercepted liquid wliq[u,s,l]max
             wsnomax, &  ! max amount of intercepted snow   wsno[u,s,l]max
             wliqmin, &  ! min amount of intercepted liquid (same name for u,s,l)
             wsnomin     ! min amount of intercepted snow (same name for u,s,l)

      real*8 tveg(npt),  &   ! temperature of veg component t[u,s,l]
             wliq(npt),  &   ! intercepted liquid amount wliq[u,s,l] (returned)
             wsno(npt),  &   ! intercepted snow amount wsno[u,s,l] (returned)
             xai(npt),   &   ! lai, sai, lai+sai for upper leaves/stems,lower veg
             pflux(npt), &   ! ht flux due to adjust of intercep precip (returned)
             rain(npt),  &   ! rainfall rate. Input: above veg, Output: below veg
             train(npt), &   ! temperature of rain. (returned)
             snow(npt),  &   ! snowfall rate. Input: above veg, output: below veg
             tsnow(npt)      ! temperature of snow (returned)

! local variables
      integer i
      real*8 rwork,      &    ! 1/dtime
             x, rwork2,  &    ! work variables
             dw               ! correction: freezing liguid or melting snow
      real*8 fint(npt),  &    ! precip fraction intercepted by unit leaf/stem area
             drip(npt),  &    ! rate of liquid drip
             blow(npt)        ! rate of snow blowoff

! ---------------------------------------------------------------------
! 
! calculate fint, the intercepted precip fraction per unit
! leaf/stem area -- note 0.5 * lai or sai (similar to irrad)
      do 50 i = 1, npt
         if (xai(i).ge.epsilon) then
            fint(i) = ( 1.-exp(-0.5*xai(i)) )/ xai(i)
         else
            fint(i) = 0.5
         endif
50    continue

! step intercepted liquid and snow amounts due to drip/blow,
! intercepted rainfall/snowfall, and min/max limits. also 
! adjust temperature of intercepted precip to current veg temp,
! storing the heat needed to do this in pflux for use in turvap
! 
! without these pfluxes, the implicit turvap calcs could not
! account for the heat flux associated with precip adjustments,
! especially changes of phase (see below), and so could not
! handle equilibrium situations such as intercepted snowfall
! being continuously melted by warm atmos fluxes, with the veg 
! temp somewhat lower than the equil atmos temp to supply heat
! that melts the incoming snow; (turvap would just change veg 
! temp to atmos equil, with little sensible heat storage...then
! final phase adjustment would return veg temp to melt point)
!
! the use of the current (ie, previous timestep's) veg temp 
! gives the best estimate of what this timestep's final temp
! will be, at least for steady conditions
      rwork = 1. / dtime
      do 100 i = 1, npt

! liquid
         drip(i) = xai(i)*wliq(i)/tdrip
         wliq(i) = wliq(i) * (1.-dtime/tdrip)
         wliq(i) = wliq(i) + dtime*rain(i)*fint(i)
         pflux(i) = rain(i)*fint(i) * (tveg(i)-train(i))*ch2o
         rain(i) = rain(i)*(1.-xai(i)*fint(i))
         x = wliq(i)
         wliq(i) = min (wliq(i), wliqmax)
         if (wliq(i).lt.wliqmin) wliq(i) = 0. 
         drip(i) = drip(i) + xai(i)*(x-wliq(i))*rwork

! snow
         blow(i) = xai(i)*wsno(i)/tblow
         wsno(i) = wsno(i) * (1.-dtime/tblow)
         wsno(i) = wsno(i) + dtime*snow(i)*fint(i)
         pflux(i) = pflux(i) + snow(i)*fint(i) * (tveg(i)-tsnow(i))*cice
         snow(i) = snow(i)*(1.-xai(i)*fint(i))
         x = wsno(i)
         wsno(i) = min (wsno(i), wsnomax)
         if (wsno(i).lt.wsnomin) wsno(i) = 0. 
         blow(i) = blow(i) + xai(i)*(x-wsno(i))*rwork
100   continue

! change phase of liquid/snow below/above melt point, and add
! required heat to pflux (see comments above). this will only
! affect the precip intercepted in this timestep, since original
! wliq, wsno must have been ge/le melt point (ensured in later
! call to cascad2/steph2o2)
      rwork2 = ch2o - cice

      do 300 i = 1, npt
! liquid below freezing
         dw = 0.
         if (tveg(i).lt.tmelt) dw = wliq(i)
         pflux(i) = pflux(i) + dw * (rwork2*(tmelt-tveg(i)) - hfus) * rwork
         wliq(i) = wliq(i) - dw
         wsno(i) = wsno(i) + dw

! snow above freezing
         dw = 0.
         if (tveg(i).gt.tmelt)  dw = wsno(i)
         pflux(i) = pflux(i) + dw * (rwork2*(tveg(i)-tmelt) + hfus) * rwork
         wsno(i) = wsno(i) - dw
         wliq(i) = wliq(i) + dw
300   continue

! adjust rainfall, snowfall below veg for interception 
! and drip, blowoff
      call mix (rain,train, rain,train, drip,tveg, vzero,vzero, npt)
      call mix (snow,tsnow, snow,tsnow, blow,tveg, vzero,vzero, npt)
      return
end subroutine steph2o
