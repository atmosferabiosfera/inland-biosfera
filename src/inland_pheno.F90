#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine pheno(kpti, kptj)
! ---------------------------------------------------------------------
     use inland_parameters
     use inland_comatm
     use inland_comsoi
     use inland_comsum
     use inland_comveg
     use inland_comcrop, only:isimagro

     implicit none
! -----------------------------------------------------------------------

! Arguments
     integer :: kpti            ! index of 1st point of little vector
                                ! in big lpt vector
     integer :: kptj            ! index of last point of little vector

! local variables
     integer :: i
     integer :: jday
     real*8 :: ddays,      & !
               ddfac,      & !
               tthreshold, & ! temperature threshold for budburst and senescence
               gthreshold, & ! gdd threshold for budburst and senescence
               avglaiu,    & ! average lai of upper canopy
               avglail,    & ! average lai of lower canopy
               sthreshold, &
               sumcrop
! -----------------------------------------------------------------------
! define 'drop days' -- number of days to affect phenology change
     ddays = 15.0
     ddfac = 1.0 / ddays

    if(isimagro .gt. 0)then
! soil temperature threshold for budburst
! soil temperature threshold is assumed to be 0 degrees C
      sthreshold = 273.16

    endif !endo of isimagro
!
! begin global grid
!

     do 100 i = kpti, kptj
!
! only do if crops are not growing
!
        if (isimagro .eq. 0) then
! ---------------------------------------------------------------------
! * * * upper canopy winter phenology * * *
! ---------------------------------------------------------------------
!
! temperature threshold for budburst and senescence
!
! temperature threshold is assumed to be 0 degrees C
! or 5 degrees warmer than the coldest monthly temperature
        tthreshold = max(dble(0.0) + 273.16, tc(i) + 5.0 + 273.16)

! gdd threshold temperature for leaf budburst
! with a growing degree threshold of 100 units
        gthreshold = 0.0 + 273.16

! determine if growing degree days are initiated
        if (a10td(i) .lt. gthreshold) then
           agddu(i)  = 0.0
        else
           agddu(i) = agddu(i) + td(i) - gthreshold
        endif

! determine leaf display
        if (a10td(i).lt.tthreshold) then
           tempu(i)  = max(dble(0.0), tempu(i) - ddfac)
        else
           tempu(i) = min(dble(1.), max(dble(0.0), agddu(i) - 100.0) / 50.0)
        endif

! ---------------------------------------------------------------------
! * * * lower canopy winter phenology * * *
! ---------------------------------------------------------------------
!
! temperature threshold for budburst and senescence
!
! temperature threshold is assumed to be 0 degrees C
        tthreshold = 0.0 + 273.16

! gdd threshold temperature for leaf budburst
! with a growing degree threshold of 150 units
        gthreshold = -5.0 + 273.16

! determine if growing degree days are initiated
        if (a10td(i).lt.gthreshold) then
           agddl(i)  = 0.0
        else
           agddl(i) = agddl(i) + td(i) - gthreshold
        endif

! determine leaf display
        if (a10td(i).lt.tthreshold) then
           templ(i)  = max(dble(0.0), templ(i) - ddfac)
        else
           templ(i) = min(dble(1.), max(dble(0.0), agddl(i) - 150.0) / 50.0)
        endif

! ---------------------------------------------------------------------
! * * * drought canopy winter phenology * * *
! ---------------------------------------------------------------------
        if (a10ancub(i).lt.0.0) dropu(i)  = max(dble(0.1), dropu(i)  - dble(ddfac))
        if (a10ancub(i).ge.0.0) dropu(i)  = min(dble(1.0), dropu(i)  + dble(ddfac))
        if (a10ancls(i).lt.0.0) dropls(i) = max(dble(0.1), dropls(i) - dble(ddfac))
        if (a10ancls(i).ge.0.0) dropls(i) = min(dble(1.0), dropls(i) + dble(ddfac))
        if (a10ancl4(i).lt.0.0) dropl4(i) = max(dble(0.1), dropl4(i) - dble(ddfac))
        if (a10ancl4(i).ge.0.0) dropl4(i) = min(dble(1.0), dropl4(i) + dble(ddfac))
        if (a10ancl3(i).lt.0.0) dropl3(i) = max(dble(0.1), dropl3(i) - dble(ddfac))
        if (a10ancl3(i).ge.0.0) dropl3(i) = min(dble(1.0), dropl3(i) + dble(ddfac))

! ---------------------------------------------------------------------
! * * * update lai and canopy fractions * * *
! ---------------------------------------------------------------------
!
! upper canopy single sided leaf area index (area-weighted)
        avglaiu = plai(i,1)             + &
                  plai(i,2) * dropu(i)  + &
                  plai(i,3)             + &
                  plai(i,4)             + &
                  plai(i,5) * tempu(i)  + &
                  plai(i,6)             + &
                  plai(i,7) * tempu(i)  + &
                  plai(i,8) * tempu(i)

! upper canopy fractions
        frac(i,1) = plai(i,1)            / max(avglaiu, epsilon)
        frac(i,2) = plai(i,2) * dropu(i) / max(avglaiu, epsilon)
        frac(i,3) = plai(i,3)            / max(avglaiu, epsilon)
        frac(i,4) = plai(i,4)            / max(avglaiu, epsilon)
        frac(i,5) = plai(i,5) * tempu(i) / max(avglaiu, epsilon)
        frac(i,6) = plai(i,6)            / max(avglaiu, epsilon)
        frac(i,7) = plai(i,7) * tempu(i) / max(avglaiu, epsilon)
        frac(i,8) = plai(i,8) * tempu(i) / max(avglaiu, epsilon)

! lower canopy single sided leaf area index (area-weighted)
        avglail = plai(i,9)                             + &
                  plai(i,10) * min(templ(i), dropls(i)) + &
                  plai(i,11) * min(templ(i), dropl4(i)) + &
                  plai(i,12) * min(templ(i), dropl3(i))

! lower canopy fractions
        frac(i,9)  = plai(i,9) / max(avglail, epsilon)
        frac(i,10) = plai(i,10) * min(templ(i),dropls(i)) / max(avglail,epsilon)
        frac(i,11) = plai(i,11) * min(templ(i),dropl4(i)) / max(avglail,epsilon)
        frac(i,12) = plai(i,12) * min(templ(i),dropl3(i)) / max(avglail,epsilon)

! calculate the canopy leaf area index using the fractional vegetation cover
        lai(i,1) = avglail / fl(i)
        lai(i,2) = avglaiu / fu(i)

! put a fix on canopy lais to avoid problems in physics
        lai(i,1) = min(lai(i,1), dble(12.0))
        lai(i,2) = min(lai(i,2), dble(12.0))

! ---------------------------------------------------------------------
! * * * update canopy height parameters * * *
! ---------------------------------------------------------------------
!
! update lower canopy height parameters
!
! note that they are based on vegetation fraction and not
! averaged over the entire gridcell
        zbot(i,1)   =  0.05
        ztop(i,1)   =  max(dble(0.25), lai(i,1) * dble(0.25))

! constrain ztop to be at least 0.5 meter lower than
! zbot for upper canopy
        ztop(i,1) = min(ztop(i,1), zbot(i,2) - 0.5)

      endif ! crop existence check 
! end of loop
100  continue
     return
end subroutine pheno
