#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine steph2o2 (tveg,wliq,wsno,fveg,xai,rliq,fvapw,fvapa, &
                     fsena,ta,cveg,npt)
! ---------------------------------------------------------------------
!
! removes evaporation from intercepted h2o, and does final
! heat-conserving adjustment for any liquid/snow below/above
! melt point, for one veg component
!
!     fsena and fvapa are passed as argument so that they can be used
!     only on the small vector.
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsatp

      implicit none
!-----------------------------------------------------------------------
! Arguments (all arguments are supplied unless otherwise noted)
      integer npt          ! number of points in little vector 
      real*8 cveg          ! specific heat of veg component ch[u,s,l] 

      real*8 tveg(npt),  & ! temperature of veg component t[u,s,l] (returned)
             wliq(npt),  & ! intercepted liquid amount wliq[u,s,l] (returned)
             wsno(npt),  & ! intercepted snow amount wsno[u,s,l] (returned)
             fveg(npt),  & ! fractional areal coverage, fu or (1-fi)*fl
             xai(npt),   & ! lai, sai, lai+sai for upper leaves/stems,lower veg
             rliq(npt),  & ! ratio of area wetted by liquid to total wetted area
             fvapw(npt), & ! wet evap h2o flx per leaf/stem area fvap[uw,s,lw]
             fvapa(npt), & ! evaporation rate
             fsena(npt), & ! sensible heat
             ta(npt)       ! air T

! local variables
      integer i

      real*8 zm,     &    ! to compute corrective fluxes
             rwork,  &    ! 1/specific heat of fusion 
             chav         ! average specific heat for veg, liw and snow
!dh(npt): correct heat flux for fluid below melt point and opposite
!dw(npt): correct water flux for fluid below melt point and opposite
      real*8 dh(npt), dw(npt)   

!-----------------------------------------------------------------------
#define STEPH2O2_COMSAT
#include "inland_comsat.h"
!-----------------------------------------------------------------------
! step intercepted h2o due to evaporation/sublimation.
! (fvapw already has been multiplied by fwet factor in turvap,
! so it is per unit leaf/stem area.)
!
! due to linear fwet factors (see comments in fwetcal) and
! the cap on suw,ssw,slw in turvap, evaporation in one timestep
! should hardly ever make wliq or wsno negative -- but if this
! happens, compensate by increasing vapor flux from atmosphere, 
! and decreasing sensib heat flux from atmos (the former is
! dangerous since it could suck moisture out of a dry atmos,
! and both are unphysical but do fix the budget) tveg in hvapf
! and hsubf should be pre-turvap-timestep values, but are not
      do 100 i = 1, npt
         wliq(i) = wliq(i) - dtime *     rliq(i)  * fvapw(i)
         wsno(i) = wsno(i) - dtime * (1.-rliq(i)) * fvapw(i)

! check to see if predicted wliq or wsno are less than zero
         if ((wliq(i).lt.0. .or. wsno(i).lt.0.) .and. &
             fveg(i)*xai(i).gt.0. )  then
!         write (*,9999) i, wliq(i), wsno(i)
!9999     format(' ***warning: wliq<0 or wsno<0 -- steph2o2 9999',
!    >           ' i, wliq, wsno:',i4, 2f12.6)
!
! calculate corrective fluxes
            zm = max (-wliq(i), dble(0.)) * fveg(i) * xai(i) / dtime
            fvapa(i) = fvapa(i) + zm
            fsena(i) = fsena(i) - zm*hvapf(tveg(i),ta(i))
            wliq(i) = max (wliq(i), dble(0.))

            zm = max (-wsno(i), dble(0.)) * fveg(i) * xai(i) / dtime
            fvapa(i) = fvapa(i) + zm
            fsena(i) = fsena(i) - zm*hsubf(tveg(i),ta(i))
            wsno(i) = max (wsno(i), dble(0.))
         endif
100   continue

! final heat-conserving correction for liquid/snow below/above
! melting point
      rwork = 1. / hfus
      do 200 i = 1, npt
         chav = cveg + ch2o*wliq(i) + cice*wsno(i)

! correct for liquid below melt point
! (nb: if tveg > tmelt or wliq = 0, nothing changes.)
         if (tveg(i).lt.tmelt .and. wliq(i).gt.0.0) then
            dh(i) = chav*(tmelt - tveg(i))
            dw(i) = min (wliq(i), max (dble(0.), dh(i)*rwork))
            wliq(i) = wliq(i) - dw(i)
            wsno(i) = wsno(i) + dw(i) 
            chav = cveg + ch2o*wliq(i) + cice*wsno(i)
            tveg(i) = tmelt - (dh(i)-hfus*dw(i))/chav
         endif

! correct for snow above melt point
!
! (nb: if tveg < tmelt or wsno = 0, nothing changes.)
         if (tveg(i).gt.tmelt .and. wsno(i).gt.0.0) then
            dh(i) = chav*(tveg(i) - tmelt)
            dw(i) = min (wsno(i), max (dble(0.), dh(i)*rwork))
            wsno(i) = wsno(i) - dw(i)
            wliq(i) = wliq(i) + dw(i)
            chav = cveg + ch2o*wliq(i) + cice*wsno(i)
            tveg(i) = tmelt + (dh(i)-hfus*dw(i))/chav
         endif
200   continue
      return
end subroutine steph2o2
