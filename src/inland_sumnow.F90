#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine sumnow(kpti, kptj)
! ---------------------------------------------------------------------
! subroutine to calculate instantaneous values of root respiration (per
! time step (dtime) and net ecosystem exchange of co2 from instantaneous
! npp calculation
! ---------------------------------------------------------------------
      use inland_parameters, only: lbeg, lend, npft, dtime, ndaypy, nsoilay,     &
                                stemae, rootae, rrootpy, rwoodpy, rgrowth
      use inland_control !, only: ipointout, iyear
      use inland_comatm
      use inland_comhyd
      use inland_comsoi
      use inland_comsno
      use inland_comsum
      use inland_comveg
      use inland_com1d
      use inland_comcrop, only: isimagro

      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      integer kpti    ! index of 1st point of little vector 
                      ! in big lpt vector
      integer kptj    ! index of last point of little vector
!------------------------------Local variables-------------------------
      integer i, k    ! loop indices
      real*8 funca    ! temperature function for aboveground biomass (stems)
      real*8 funcb    ! temperature function for belowground biomass (roots)
      real*8 roottemp ! average root temperature for all roots
!      real*8 rgrowth  ! growth respiration coefficient (fraction)
      real*8 rgrowthc ! growth respiration coefficient (fraction) for crops
      real*8 rroot    ! maintenance respiration coefficient for root (/s)
      real*8 rwood    ! maintenance respiration coefficient for wood (/s)
      real*8 stemtemp ! stem temperature
      real*8 zweight  ! 10-day time averaging factor
      real*8 workday
      real*8 smask    ! 1 - fi
      real*8 xminlai  ! minimum LAI used as seed when climate becomes favorable
      real*8 rwork    ! kgC m-2 yr-1 to mol CO2 m-2 s-1
      real*8 zweight11! 11-day time averaging factor
      real*8 zweight3 ! 3-day time averaging factor
      real*8 zweight5 ! 5-day time averaging factor

      real*8 cbiolmin(lbeg:lend,npft) ! minimum leaf biomass sed as seed.
      real*8 nppdummy(lbeg:lend,npft) ! temporary variable to calculate npp
! ---------------------------------------------------------------------
! * * * define working variables * * *
! ---------------------------------------------------------------------

! maintenance respiration coefficients (per second)
!
! initially, we pick values for respiration coefficients that
! defined in units of  / year
!
!   rwood ~ 0.0125 
!   rroot ~ 1.2500
!
! however, we convert the unitsconvert to have resulting respiration
! fluxes in units of mol-CO2 / m**2 / second
!
! this requires we convert the time unit to seconds and add an additional
! factor to convert biomass units from kilograms to moles
      rwood   = rwoodpy / (ndaypy * 86400.0) * (1000.0 / 12.0)
      rroot   = rrootpy / (ndaypy * 86400.0) * (1000.0 / 12.0)

! growth respiration coefficient (fraction)
     if(isimagro .gt. 0) then
      rgrowth = 0.30
      rgrowthc = 0.242  ! crop growth respiration for sugarcane from Singles,2002
     endif

! number of timesteps per day (for use in gdd0this)
      workday = 86400./dtime

! 10-day time averaging factor
      zweight = exp(-1. / (10.0 * workday))
  if(isimagro .gt. 0) then
!
! 11-day time averaging factor
!
      zweight11 = exp(-1. / (11.0 * workday))
!
! 3-day time averaging factor
!
      zweight3 = exp(-1. / (3.0 * workday))
!
! 5-day time averaging factor
!
      zweight5 = exp(-1. / (5.0 * workday)) 
!
  endif

!Commented by Pousa, it doesn' necessary
! factors for instantaneous NEE 
!      rwork = 1./(ndaypy * 86400.0 * 12.e-3 )     
!      xminlai =  0.010

! begin global grid
      do 100 i = kpti, kptj

! calculate instantaneous carbon flux parameters, including
! npp (net primary production) and nee (net ecosystem exchange)
!
! in this routine, all of the fluxes are calculated in the units
! of mol-CO2 / m**2 / sec
!
! ---------------------------------------------------------------------
! * * * calculate instantaneous GPP * * *
! ---------------------------------------------------------------------
!
! snow masking for lower canopy vegetation
         smask = 1.0 - fi(i)

! note that the following plants types follow different physiological paths
!
!   - broadleaf trees   :  types 1, 2, 3, 5, 7, 8 
!   - conifer   trees   :  types 4, 6
!   - shrubs            :  types 9, 10
!   - c4 grasses        :  type 11
!   - c3 grasses        :  type 12
!   - c3 crops          :  soybean, wheat (type 13, 15) 
!   - c4 crops          :  maize and sugarcane (type 14, 16)  
!
! note that plant type 8 is actually a deciduous conifer (e.g., Larix), but
! we are assuming that it's physiological behavior is like a broadleaf tree
!
! nppdummy is canopy npp before accounting for stem & root respiration
! Navin Sept 02
         nppdummy(i,1)  = frac(i,1)  * ancub(i) * lai(i,2) * fu(i)
         nppdummy(i,2)  = frac(i,2)  * ancub(i) * lai(i,2) * fu(i)
         nppdummy(i,3)  = frac(i,3)  * ancub(i) * lai(i,2) * fu(i)
         nppdummy(i,4)  = frac(i,4)  * ancuc(i) * lai(i,2) * fu(i)
         nppdummy(i,5)  = frac(i,5)  * ancub(i) * lai(i,2) * fu(i)
         nppdummy(i,6)  = frac(i,6)  * ancuc(i) * lai(i,2) * fu(i)
         nppdummy(i,7)  = frac(i,7)  * ancub(i) * lai(i,2) * fu(i)
         nppdummy(i,8)  = frac(i,8)  * ancub(i) * lai(i,2) * fu(i)
         nppdummy(i,9)  = frac(i,9)  * ancls(i) * lai(i,1) * fl(i) * smask 
         nppdummy(i,10) = frac(i,10) * ancls(i) * lai(i,1) * fl(i) * smask
         nppdummy(i,11) = frac(i,11) * ancl4(i) * lai(i,1) * fl(i) * smask
         nppdummy(i,12) = frac(i,12) * ancl3(i) * lai(i,1) * fl(i) * smask

#ifndef SINGLE_POINT_MODEL
       if(isimagro .gt. 0)then
         nppdummy(i,13) = frac(i,13) * ancc3(i) * lai(i,1) * fl(i) * smask
         nppdummy(i,14) = frac(i,14) * ancc4(i) * lai(i,1) * fl(i) * smask
         nppdummy(i,15) = frac(i,15) * ancc3(i) * lai(i,1) * fl(i) * smask
         nppdummy(i,16) = frac(i,16) * ancc4(i) * lai(i,1) * fl(i) * smask
       endif
#endif /* SINGLE_POINT_MODEL */

! Navin's correction to compute npp using tgpp via agXXX
! agXXX should be used 
         tgpp(i,1)  = frac(i,1)  * agcub(i) * lai(i,2) * fu(i)
         tgpp(i,2)  = frac(i,2)  * agcub(i) * lai(i,2) * fu(i)
         tgpp(i,3)  = frac(i,3)  * agcub(i) * lai(i,2) * fu(i)
         tgpp(i,4)  = frac(i,4)  * agcuc(i) * lai(i,2) * fu(i)
         tgpp(i,5)  = frac(i,5)  * agcub(i) * lai(i,2) * fu(i)
         tgpp(i,6)  = frac(i,6)  * agcuc(i) * lai(i,2) * fu(i)
         tgpp(i,7)  = frac(i,7)  * agcub(i) * lai(i,2) * fu(i)
         tgpp(i,8)  = frac(i,8)  * agcub(i) * lai(i,2) * fu(i)
         tgpp(i,9)  = frac(i,9)  * agcls(i) * lai(i,1) * fl(i) * smask 
         tgpp(i,10) = frac(i,10) * agcls(i) * lai(i,1) * fl(i) * smask
         tgpp(i,11) = frac(i,11) * agcl4(i) * lai(i,1) * fl(i) * smask
         tgpp(i,12) = frac(i,12) * agcl3(i) * lai(i,1) * fl(i) * smask
       if(isimagro .gt. 0)then
         tgpp(i,13) = frac(i,13) * agcc3(i) * lai(i,1) * fl(i) * smask
         tgpp(i,14) = frac(i,14) * agcc4(i) * lai(i,1) * fl(i) * smask
         tgpp(i,15) = frac(i,15) * agcc3(i) * lai(i,1) * fl(i) * smask
         tgpp(i,16) = frac(i,16) * agcc4(i) * lai(i,1) * fl(i) * smask
       endif

! calculate total gridcell gpp
         tgpptot(i) = 0.0
         do 110 k = 1, npft
            tgpptot(i) = tgpptot(i) + tgpp(i,k)
110      continue

! ---------------------------------------------------------------------
! * * * calculate temperature functions for respiration * * *
! ---------------------------------------------------------------------
!
! calculate the stem temperature
         stemtemp = ts(i)

! calculate average root temperature (average of all roots)
         roottemp = 0.0
         do 120 k = 1, nsoilay
            roottemp = roottemp + tsoi(i,k) * 0.5 * (froot(k,1) + froot(k,2))
120      continue

! calculate respiration terms on a 15 degree base
! following respiration parameterization of Lloyd and Taylor
         funca = exp(stemae * (1. / 288.16 - 1. / stemtemp))
         funcb = exp(rootae * (1. / 288.16 - 1. / roottemp))

! ---------------------------------------------------------------------
! * * * calculate instantaneous NPP * * *
! ---------------------------------------------------------------------
!
! the basic equation for npp is
!
!   npp = (1 - growth respiration term) * (gpp - maintenance respiration terms)
!
! here the respiration terms are simulated as
!
!   growth respiration = rgrowth * (gpp - maintenance respiration terms)
!
! where
!
!   rgrowth is the construction cost of new tissues
!
! and
!
!   root respiration = rroot * cbior(i,k) * funcb
!   wood respiration = rwood * cbiow(i,k) * funca * sapwood fraction
!
! where
! 
!   funca = temperature function for aboveground biomass (stems)
!   funcb = temperature function for belowground biomass (roots)
!
! note that we assume the sapwood fraction for shrubs is 1.0
!
! also note that we apply growth respiration, (1 - rgrowth), 
! throughout the year; this may cause problems when comparing
! these npp values with flux tower measurements
!
! also note that we need to convert the mass units of wood and
! root biomass from kilograms of carbon to moles of carbon
! to maintain consistent units (done in rwood, rroot)
!
! finally, note that growth respiration is only applied to 
! positive carbon gains (i.e., when gpp-rmaint is positive)
         tnpp(i,1)  = nppdummy(i,1) - rwood * cbiow(i,1) * sapfrac(i) * funca -&
                      rroot * cbior(i,1) * funcb
         tnpp(i,2)  = nppdummy(i,2) - rwood * cbiow(i,2) * sapfrac(i) * &
                     funca - rroot * cbior(i,2) * funcb
         tnpp(i,3)  = nppdummy(i,3) - rwood * cbiow(i,3) * sapfrac(i) * &
                      funca - rroot * cbior(i,3) * funcb
         tnpp(i,4)  = nppdummy(i,4) - rwood * cbiow(i,4) * sapfrac(i) * &
                      funca - rroot * cbior(i,4) * funcb
         tnpp(i,5)  = nppdummy(i,5) - rwood * cbiow(i,5) * sapfrac(i) * &
                      funca - rroot * cbior(i,5) * funcb
         tnpp(i,6)  = nppdummy(i,6) - rwood * cbiow(i,6) * sapfrac(i) * &
                      funca - rroot * cbior(i,6) * funcb
         tnpp(i,7)  = nppdummy(i,7) - rwood * cbiow(i,7) * sapfrac(i) * &
                      funca -  rroot * cbior(i,7) * funcb
         tnpp(i,8)  = nppdummy(i,8) - rwood * cbiow(i,8) * sapfrac(i) * &
                      funca - rroot * cbior(i,8) * funcb
         tnpp(i,9)  = nppdummy(i,9) - rwood * cbiow(i,9) * funca - rroot * &
                      cbior(i,9) * funcb
         tnpp(i,10) = nppdummy(i,10) - rwood * cbiow(i,10) * funca - rroot * &
                      cbior(i,10) * funcb
         tnpp(i,11) = nppdummy(i,11) - rroot * cbior(i,11) * funcb
         tnpp(i,12) = nppdummy(i,12) - rroot * cbior(i,12) * funcb
#ifndef SINGLE_POINT_MODEL
       if(isimagro .gt. 0)then
         tnpp(i,13) = nppdummy(i,13) - rroot * cbior(i,13) * funcb
         tnpp(i,14) = nppdummy(i,14) - rroot * cbior(i,14) * funcb
         tnpp(i,15) = nppdummy(i,15) - rroot * cbior(i,15) * funcb
         tnpp(i,16) = nppdummy(i,16) - rroot * cbior(i,16) * funcb

!FIXME Variable below is not being used anywere else
!      	  rootr(i) = rroot * cbior(i,16) * funcb   !sugarcane 
       endif
#endif /* SINGLE_POINT_MODEL */
! Apply existence array from previous year to have coherent instanteneous nee and yearly nee.
!         do k = 1, npft
!            tnpp(i,k) = tnpp(i,k) * exist(i,k)
!         end do
!
         if(isimagro .eq. 0) then
          do k = 1, npft
            tnpp(i,k) = tnpp(i,k) * exist(i,k)
          end do
         endif

! apply growth respiration and calculate total gridcell npp
!
        tnpptot(i) = 0.0
!
        do 130 k = 1, npft
!
! CK non-crop types
!
          if (tnpp(i,k).gt.0.0 .and. k .le. 12) &
             tnpp(i,k) = tnpp(i,k)  * (1.0 - rgrowth)
! CK : crop types
!
      if(isimagro .gt. 0)then
          if (tnpp(i,k).gt.0.0 .and. k .gt. 12) &
             tnpp(i,k) = tnpp(i,k)  * (1.0 - rgrowthc)
      endif
          tnpptot(i) = tnpptot(i) + tnpp(i,k)
 130    continue
!
! ---------------------------------------------------------------------
! * * * calculate total fine root respiration * * *
! ---------------------------------------------------------------------
         tco2root(i) = 0.0
         do 140 k = 1, npft
            tco2root(i) = tco2root(i) + rroot * cbior(i,k) * funcb
140      continue

! ---------------------------------------------------------------------
! * * * calculate instantaneous NEE * * *
! ---------------------------------------------------------------------
!
! microbial respiration is calculated in biogeochem.f
!
! Units are (mol CO2 / m^2 / timestep)
         tneetot(i) = tnpptot(i) - tco2mic(i)

! ---------------------------------------------------------------------
! * * * update 10-day running-mean parameters * * *
! ---------------------------------------------------------------------
!
! 10-day daily air temperature
         a10td(i)    = zweight * a10td(i)    + (1. - zweight) * td(i)

     if(isimagro .gt. 0) then
        a10ts(i) = zweight  * a10ts(i) + (1. - zweight)  * tsoi(i,1)
        a5td(i)  = zweight5 * a5td(i)  + (1. - zweight5) * td(i)
!
! 3-day daily minimum air temperature
!
        a3tdmin(i) = zweight3 * a3tdmin(i) + (1.-zweight3) * tmin(i)
!
! 11-day running mean daily soil temperature
! Found after an estimate of soil temperature at 10cm from top two
! IBIS soil layers (0-10 and 10-25cm) is determined
!
        tsoiavg(i) = 0.6 * tsoi(i,1) + 0.4 * tsoi(i,2)
        a11soiltd(i) = zweight11 * a11soiltd(i) + (1.-zweight11) * tsoiavg(i)
!
     endif
! 10-day canopy photosynthesis rates
         a10ancub(i) = zweight * a10ancub(i) + (1. - zweight) * ancub(i)
         ! if this is not needed by w/r restart, then it is not needed at all!
         !a10ancuc(i) = zweight * a10ancuc(i) + (1. - zweight) * ancuc(i)
         a10ancls(i) = zweight * a10ancls(i) + (1. - zweight) * ancls(i)
         a10ancl3(i) = zweight * a10ancl3(i) + (1. - zweight) * ancl3(i)
         a10ancl4(i) = zweight * a10ancl4(i) + (1. - zweight) * ancl4(i)
     if(isimagro .gt. 0) then
         a10ancc3(i) = zweight * a10ancc3(i) + (1. - zweight) * ancc3(i)
         a10ancc4(i) = zweight * a10ancc4(i) + (1. - zweight) * ancc4(i)

! 10-day minimimum daily temperature average -- used to help determine
! planting dates of crops.  If both the daily average temperature
! (10-day mean) > 8 C and the 10-day minimum temperature average
! is > 0 C, then soybean and maize crops are allowed to be planted. 
!
        a5tmin(i)  = zweight5 * a5tmin(i) + (1. - zweight5) * tmin(i)
        a10tmin(i) = zweight  * a10tmin(i)+ (1. - zweight)  * tmin(i)
     endif
100   continue
      if (ipointout .eq. 1) then
         i = 1001
         if (i .ge. kpti .and. i .le. kptj) then
            write(121,'(13f15.3)') tnpptot(i)*1e9,(tnpp(i,k)*1e9,k=1,12)
!           call flush(121)
         end if
         i = 662
         if (i .ge. kpti .and. i .le. kptj) then
            write(1121,'(13f15.3)') tnpptot(i)*1e9,(tnpp(i,k)*1e9,k=1,12)
!           call flush(1121)
         end if
         i = 651
         if (i .ge. kpti .and. i .le. kptj) then
            write(2121,'(13f15.3)') tnpptot(i)*1e9,(tnpp(i,k)*1e9,k=1,12)
!           call flush(2121)
         end if
      end if
      return
end subroutine sumnow
