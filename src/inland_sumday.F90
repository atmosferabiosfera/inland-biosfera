#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine sumday (mcsec, loopi, kpti, kptj)
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control
      use inland_comatm
      use inland_comhyd
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_com1d
      use inland_comfire
      use inland_comcrop, only:isimagro
      use inland_comhour, only:imetyear

      implicit none
!-----------------------------------------------------------------------
! input-output variables
      integer mcsec           ! current second in the day (passed in, it's time)
      integer loopi           ! index of little vector in big vector
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer i, k, np, ndtimes
      real*8 rwork, rwork2, rwork3, rwork4, & ! work time variable
             tconst,   & ! constant for Lloyd and Taylor (1994) function
             bconst,   & ! base temperature used for carbon decomposition
             btemp,    & ! maximum value of decomposition factor
             rdepth,   & ! total depth of the 4 1st soil layers
             rdepth2,  & ! total depth of the 2 1st soil layers
             snodpth,  & ! total snow depth
             soiltemp, & ! average soil temp for 2 1st layers
             soilmois, & ! average soil moisture (fraction of porosity) for 2 1st layers
             soilice,  & ! average soil ice for 2 1st layers
             soitempc, & ! average soil temp over 6 layers
             soimoisc, & ! average soil moisture over 6 layers
             factor,   & ! temperature decompositn factor for liter/soil carbon
             wfps,     & ! water filled pore space
             moist,    & ! moisture effect on decomposition
             awc,      & ! available water content (fraction)
             soilmois2,& ! average soil ice for 2 1st layers
             precipfac   !

! ---------------------------------------------------------------------
! * * * update counters and working variables * * *
! ---------------------------------------------------------------------
!
! reset sumday if the first timestep of the day
! FIXME: dangerous test. It should be mcsec .ge. dtime!
      if (mcsec .eq. 0) then
         ndtime(loopi) = 0
         do i = kpti, kptj
           if(isimagro .eq. 0)then
            gdd0this(i) = gdd0this(i) + max(dble(0.), (td(i) - 273.16))
            gdd5this(i) = gdd5this(i) + max(dble(0.), (td(i) - 278.16))
           endif
         end do
         adburnfrac(:) = 0
      end if

! accumulate daily output (at this point for soil decomposition)
      ndtime(loopi) = ndtime(loopi) + 1

! put into local integer variable
      ndtimes = ndtime(loopi)

! working variables
      rwork  = 1. / float(ndtimes)
      rwork2 = 86400.
      rwork3 = 86400. * 12.e-3
      rwork4 = 86400. * 14.e-3

! constants used in temperature function for c decomposition
! (arrhenius function constant)
      tconst  = 344.00  ! constant for Lloyd and Taylor (1994) function
      btemp   = 288.16  ! base temperature used for carbon decomposition
      bconst  = 10.0    ! maximum value of decomposition factor

! soil weighting factors
      rdepth  = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
      if(isimagro .eq. 0)then
         rdepth2 = 1. / (hsoi(1) + hsoi(2))
      else
         rdepth2 = 1. / (hsoi(1) + hsoi(2)+ hsoi(3))
      endif
!
! begin global grid
      do 100 i = kpti, kptj
      if(isimagro .gt. 0) then
!
! calculation of daily npp values for crops
!
        adnpp(i,13) = ((ndtimes-1) * adnpp(i,13) + tnpp(i,13) * rwork3) * rwork
        adnpp(i,14) = ((ndtimes-1) * adnpp(i,14) + tnpp(i,14) * rwork3) * rwork
        adnpp(i,15) = ((ndtimes-1) * adnpp(i,15) + tnpp(i,15) * rwork3) * rwork
        adnpp(i,16) = ((ndtimes-1) * adnpp(i,16) + tnpp(i,16) * rwork3) * rwork
        adnpp(i,17) = ((ndtimes-1) * adnpp(i,17) + tnpp(i,17) * rwork3) * rwork
!
      endif ! end the isimagro
! ---------------------------------------------------------------------
! * * * daily water budget terms * * *
! ---------------------------------------------------------------------
         adrain(i)    = ((ndtimes-1) * adrain(i) + raina(i) * 86400.) * rwork
         adsnow(i)    = ((ndtimes-1) * adsnow(i) + snowa(i) * 86400.) * rwork
         adaet(i)     = ((ndtimes-1) * adaet(i) - fvapa(i) * 86400.) * rwork
         adtrunoff(i) = ((ndtimes-1) * adtrunoff(i) + (grunof(i) + gdrain(i)) *&
                        86400.) * rwork
         adsrunoff(i)  = ((ndtimes-1) * adsrunoff(i) + grunof(i) * 86400.) *   &
                         rwork
         addrainage(i) = ((ndtimes-1) * addrainage(i) + gdrain(i) * 86400.) *  &
                         rwork

      if(isimagro .gt. 0) then
        adtrans(i)   = ((ndtimes-1) * adtrans(i) + (gtransl(i) + gtransu(i)) * 86400.) * rwork
        adevap(i)    = adaet(i) - adtrans(i)
        adtratio(i)  = max(0.0, min(1.0, adtrans(i) / adaet(i)))

        if(istep.eq.24)write(226,*),iyear,jday, adrain(i)

      endif


! ---------------------------------------------------------------------
! * * * daily atmospheric terms * * *
! ---------------------------------------------------------------------
!
      if(isimagro .eq. 0) then
! Different from off-line IBIS where td comes from climatology
! Daily mean temperature used for phenology based on 2-m screen
! temperature instead of 1st atmospheric level (~ 70 m)
!        td(i)      = ((ndtimes-1) * td(i) + ta(i)) * rwork
         td(i)   = ((ndtimes-1) * td(i) + ts2(i)) * rwork
      else

!        adrh(i) = ((ndtimes-1) * adrh(i) + rh(i)) * rwork
!        adud(i) = ((ndtimes-1) * adud(i) + ud(i)) * rwork
      endif
!
! ---------------------------------------------------------------------
! * * * daily fire parameters * * *
! ---------------------------------------------------------------------
if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
         adpbio(i) = ((ndtimes-1) * adpbio(i) + pbio(i)) * rwork
         adpmoi(i) = ((ndtimes-1) * adpmoi(i) + pmoi(i)) * rwork
         adpign(i) = ((ndtimes-1) * adpign(i) + pign(i)) * rwork
         adpfire(i) = ((ndtimes-1) * adpfire(i) + pfire(i)) * rwork
         adsrate(i) = ((ndtimes-1) * adsrate(i) + srate(i)) * rwork
         adabday(i) = ((ndtimes-1) * adabday(i) + abday(i)) * rwork
         adburnfrac(i) = ((ndtimes-1) * adburnfrac(i) + burnfrac(i)) * rwork

         do np = 1, npft
            adburnpft(i,np) = ((ndtimes-1) * adburnpft(i,np) + burnpft(i,np)) * rwork
         end do
end if
! ---------------------------------------------------------------------
! * * * daily snow parameters * * *
! ---------------------------------------------------------------------
         snodpth = hsno(i,1) + hsno(i,2) + hsno(i,3)
         adsnod(i) = ((ndtimes-1) * adsnod(i) + snodpth) * rwork
         adsnof(i) = ((ndtimes-1) * adsnof(i) + fi(i))   * rwork

! ---------------------------------------------------------------------
! * * * soil parameters * * *
! ---------------------------------------------------------------------
!
! initialize average soil parameters
         soiltemp = 0.0
         soilmois = 0.0
         soilmois2 = 0.0
         soilice  = 0.0
         soitempc = 0.0
         soimoisc = 0.0

         awc      = 0.0

! averages for first 2 layers of soil
         do 110 k = 1, 3
            soiltemp =  soiltemp + tsoi(i,k)  * hsoi(k)
            soilmois =  soilmois + wsoi(i,k)  * hsoi(k)
            soilice  =  soilice  + wisoi(i,k) * hsoi(k)
110      continue

! weighting on just thickness of each layer
         soilmois = soilmois * rdepth2
         soilice  = soilice  * rdepth2
         soiltemp = soiltemp * rdepth2

      if(isimagro .gt. 0) then
        do  k = 4, 6
          soilmois2 =  soilmois2 + wsoi(i,k)  * hsoi(k)
	enddo

          soilmois2 = soilmois2 / (hsoi(4)+hsoi(5)+hsoi(6))
      endif !isimagro

!
! calculate average root temperature, soil temperature and moisture and
! ice content based on rooting profiles (weighted) from jackson et al
! 1996
!
! these soil moisture and temperatures are used in biogeochem.f
! we assume that the rooting profiles approximate
! where carbon resides in the soil
         do 120 k = 1, nsoilay
          soitempc = soitempc + tsoi(i,k) * 0.5 * (froot(k,1) + froot(k,2))
          soimoisc = soimoisc + wsoi(i,k) * 0.5 * (froot(k,1) + froot(k,2))
 120     continue

! calculate daily average soil moisture and soil ice
! using thickness of each layer as weighting function
         adwsoi(i)  = ((ndtimes-1) * adwsoi(i)  + soilmois) * rwork
         adwsoi2(i)  = ((ndtimes-1) * adwsoi2(i)  + soilmois2) * rwork
         adtsoi(i)  = ((ndtimes-1) * adtsoi(i)  + soiltemp) * rwork
         adwisoi(i) = ((ndtimes-1) * adwisoi(i) + soilice)  * rwork
!
     if(isimagro .gt. 0) then
! calculate daily average soil moisture, ice, and temperature for
! each soil layer
!
! also calculate daily uptake of water by plant (mm/day)
! calculate daily average nitrate concentration in solution (mg/liter)
!
        do 130 k = 1, nsoilay

          adwsoilay(i,k)  = ((ndtimes-1)  * adwsoilay(i,k) + wsoi(i,k))  * rwork
          adwisoilay(i,k) = ((ndtimes-1)  * adwisoilay(i,k) + wisoi(i,k)) * rwork
!          adtsoilay(i,k)  = ((ndtimes-1)  * adtsoilay(i,k) + tsoi(i,k))  * rwork
!
!          adupsoil(i,k)   = ((ndtimes-1)  * adupsoil(i,k) + upsoil(i,k) * 86400.) * rwork
!
          adcsoln(i,k)    = ((ndtimes-1)  * adcsoln(i,k) + csoln(i,k)) * rwork
!
 130   continue
     endif ! check for crop existence
!
! calculate daily average for soil temp/moisture of top layer
         adtlaysoi(i) = ((ndtimes-1) * adtlaysoi(i) + tsoi(i,1)) * rwork
         adwlaysoi(i) = ((ndtimes-1) * adwlaysoi(i) + wsoi(i,1)) * rwork

!        if (istep .eq. 24) then

!          if((iyear.eq.2009.and.jday.ge.353).or.(iyear.eq.2010.and.jday.le.115)) then
          if(imetyear .ne. 9999) then

            write(230,337)iyear, jday, adrain(i) ,adaet(i),adevap(i),adtrans(i)          &
                         ,adwsoi(i)*poros(1,1)*20.*10., adwsoi2(i)*poros(1,1)*30.*10.    &
                         ,(adwsoi(i)*0.4 + adwsoi2(i)*0.6)*poros(1,1)*50.*10.            &
                         ,adtrunoff(i),adsrunoff(i),addrainage(i)

!          endif

	endif

 337    format (2(i4,1x),11(1x,f6.2))

! calculate separate variables to keep track of weighting using
! rooting profile information
!
! note that these variables are only used for diagnostic purposes
! and that they are not needed in the biogeochemistry code
!
        if(isimagro .gt. 0) adwsoic(i)  = ((ndtimes-1) * adwsoic(i) + soimoisc) * rwork
        adtsoic(i)  = ((ndtimes-1) * adtsoic(i) + soitempc) * rwork
!
! ---------------------------------------------------------------------
! * * * calculate daily soil co2 fluxes * * *
! ---------------------------------------------------------------------
!
! increment daily total co2 respiration from microbes
! tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
         adco2mic(i) = ((ndtimes-1) * adco2mic(i) + &
                       tco2mic(i) * rwork3) * rwork

! increment daily total co2 respiration from fine roots
! tco2root is instantaneous value of co2 flux calculated in stats.f
         adco2root(i) = ((ndtimes-1) * adco2root(i) + &
                        tco2root(i) * rwork3) * rwork

! calculate daily total co2 respiration from soil
!

         adco2soi(i)  = adco2root(i) + adco2mic(i)

! calculate daily ratio of total root to total co2 respiration
         if (adco2soi(i).gt.0.0) then
            adco2ratio(i) = adco2root(i) / adco2soi(i)
         else
            adco2ratio(i) = -999.99
         endif

! ---------------------------------------------------------------------
! * * * calculate daily litter decomposition parameters * * *
! ---------------------------------------------------------------------
!
! calculate litter carbon decomposition factors
! using soil temp, moisture and ice for top soil layer
!
! calculation of soil biogeochemistry decomposition factors
! based on moisture and temperature affects on microbial
! biomass dynamics
!
! moisture function based on water-filled pore space (wfps)
! williams et al., 1992 and friend et al., 1997 used in the
! hybrid 4.0 model; this is based on linn and doran, 1984
!
! temperature functions are derived from arrhenius function
! found in lloyd and taylor, 1994 with a 15 c base
!
! calculate temperature decomposition factor
         if (tsoi(i,1) .gt. 237.13) then
            factor = min (exp(tconst * ((1. / (btemp - 227.13)) - &
                     (1. / (tsoi(i,1) - 227.13)))), bconst)
         else
            factor = exp(tconst * ((1. / (btemp - 227.13)) - (1. / &
                     (237.13-227.13))))
         endif

! calculate water-filled pore space (in percent)
!
! wsoi is relative to pore space not occupied by ice and water
! thus must include the ice fraction in the calculation
         wfps = (1.0 - wisoi(i,1)) * wsoi(i,1) * 100.0

! calculate moisture decomposition factor
         if (wfps .ge. 60.0) then
            moist = 0.000371 * (wfps**2) - (0.0748 * wfps) + 4.13
         else
            moist = exp((wfps - 60.0)**2 / (-800.0))
         endif

! calculate combined temperature / moisture decomposition factor
         factor = max (dble(0.001), min (bconst, factor * moist))

! calculate daily average litter decomposition factor
         decompl(i) = ((ndtimes-1) * decompl(i) + factor) * rwork

! ---------------------------------------------------------------------
! * * * calculate daily soil carbon decomposition parameters * * *
! ---------------------------------------------------------------------
!
! calculate soil carbon decomposition factors
! using soil temp, moisture and ice weighted by rooting profile scheme
!
! calculation of soil biogeochemistry decomposition factors
! based on moisture and temperature affects on microbial
! biomass dynamics
!
! moisture function based on water-filled pore space (wfps)
! williams et al., 1992 and friend et al., 1997 used in the
! hybrid 4.0 model; this is based on linn and doran, 1984
!
! temperature functions are derived from arrhenius function
! found in lloyd and taylor, 1994 with a 15 c base
!
! calculate temperature decomposition factor
         if (soiltemp .gt. 237.13) then
            factor = min (exp(tconst * ((1. / (btemp - 227.13)) - &
                     (1. / (soiltemp - 227.13)))), bconst)
         else
            factor = exp(tconst * ((1. / (btemp - 227.13)) - (1. / &
                     (237.13-227.13))))
         endif

! calculate water-filled pore space (in percent)
!
! wsoi is relative to pore space not occupied by ice and water
! thus must include the ice fraction in the calculation
         wfps = (1. - soilice) * soilmois * 100.0

! calculate moisture decomposition factor
         if (wfps .ge. 60.0) then
            moist = 0.000371 * (wfps**2) - (0.0748 * wfps) + 4.13
         else
            moist = exp((wfps - 60.0)**2 / (-800.0))
         endif

! calculate combined temperature / moisture decomposition factor
         factor = max (dble(0.001), min (bconst, factor * moist))

! calculate daily average soil decomposition factor
         decomps(i) = ((ndtimes-1) * decomps(i) + factor) * rwork

! ---------------------------------------------------------------------
! * * * calculate other daily biogeochemical parameters * * *
! ---------------------------------------------------------------------
!
! increment daily total of net nitrogen mineralization
! value for tnmin is calculated in biogeochem.f
         adnmintot(i) = ((ndtimes-1) * adnmintot(i) + tnmin(i) * rwork4) * rwork
100   continue

      return
end subroutine sumday
