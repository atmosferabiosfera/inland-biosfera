#include "inland_config.h"
! ---------------------------------------------------------------------------------------------
! This subroutine calculates:
! 1. fire probability due to natural (lightning) and human causes
! constrained by biomass availability and soil moisture
! 2. spread rate of fire and total burned area
!
! This fire model is based on CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) V1.0
! disturbance subroutine
!-------------------------------------------------------------------------------
subroutine firectem(kpti,kptj)
!-----------------------------------------------------------------------------------------

! Global variables

      use inland_parameters, only: npft, nsoilay, garea
      use inland_comveg,   only: cbiow, cbiol, totalit, frac, fu, fl
      use inland_comsoi,   only: wsoi, swilt, sfield, wisoi
      use inland_comfire
      use inland_comatm,   only: ua
      use inland_comwork,  only: seedvec

      implicit none

! Local variables

      integer kpti, kptj, & ! initial and final number of points
              i, j          ! counters

      real Pligh            ! fire probability dependence on lightning (fraction)
      real ran_uni          ! random number generator function

      real*8 bag(kptj),   &    ! sum of aboveground biomass from leaf, stem and litter 
             fsum(kptj),  &    ! total vegetated fraction (stem + leaf)
             Ifire(kptj), &    ! indicate ignition of fire if equals to 1
             betas(kptj), &    ! beta moisture for calculating fire spread rate
             hs(kptj),    &    ! soil moisture function used for fire spread rate
             wind(kptj),  &    ! wind speed in km/h
             ws(kptj),    &    ! wind function for fire spread rate
             areax(kptj), &    ! multiplier of abday to estimate average area burned
	     sfrac(kptj), &    ! sum of upper and lower canopy fractions
	     cfrac(kptj), &    ! coefficient for total vegetated area calculation
	     sburnpft(kptj), & ! sum of burned fraction in each PFT
	     burnlow(kptj), &  ! burned fraction in lower canopy
	     burnup(kptj)      ! burned fraction in upper canopy
	       

      real*8 bsl(kptj,npft),   & ! sum of biomass from stem and leaf for PFT and land point
             betad(kptj,npft), & ! dryness term for soil layers
             wsoil(kptj,npft)    ! sum of liquid and frozen soil moisture content

      real, parameter :: PI=3.1415926535898d0

! Variables initialization

      bag(:) = 0.     
      fsum(:) = 0.    
      betas(:) = 0.   
      hs(:) = 0.      
      wind(:) = 0.    
      ws(:) = 0.
      vegfrac(:) = 0.
      sfrac(:) = 0.
      cfrac(:) = 0.
      sburnpft(:) = 0.
      burnlow(:) = 0.
      burnup(:) = 0.

      bsl(:,:) = 0.   
      betad(:,:) = 0. 
      wsoil(:,:) = 0. 

!--------------------------------------------------------------------------------------------------!

!-------------------------------------------------------
!   Probability of fire P = Pbio * Pmoi * Pign
!
! 1. Dependence on total biomass, Pbio
! 2. Dependence on soil moisture, Pmoi
! 3. Dependence on ignition sources, Pign
! 4. Probability of fire
!-------------------------------------------------------

!-------------------------------------
! 1. Dependence on total biomass, Pbio

! Sum of biomass from stem (cbiow) and leaf (cbiol) for each PFT and land point

      do 210 j=1,npft
         do 200 i=kpti,kptj
            bsl(i,j) = cbiow(i,j) + cbiol(i,j)
200	 continue
210   continue

! Sum of biomass from stem and leaf over the vegetated fraction

      do 230 j=1,npft
         do 231 i=kpti,kptj
            bag(i) = bag(i) + bsl(i,j)*frac(i,j)
231	 continue
230   continue

! Averaged biomass over vegetated fraction

      do 232 i=kpti,kptj
         fsum(i) = frac(i,1) + frac(i,2) + frac(i,3) + frac(i,4) + frac(i,5) + &
                   frac(i,6) + frac(i,7) + frac(i,8) + frac(i,9) + frac(i,10) + &
                   frac(i,11) + frac(i,12)
         if (fsum(i) .gt. 0.0) then
            bag(i) = bag(i) / fsum(i)
         else
            bag(i) = 0.0
         endif
232   continue

! Sum of biomass from stem, leaf and litter (totalit)

      do 240 i=kpti,kptj
         bag(i) = bag(i) + totalit(i)
240   continue

!   Biomass fire probability term, Pbio
! = 0, bag <= blow
! = 1 , bag >= bup
! = (bag-blow)/(bup-blow), blow < bag < bup

      do 250 i=kpti,kptj
         if (bag(i) .ge. bup) then
            Pbio(i) = 1.0
         else if (bag(i) .lt. bup .and. bag(i) .gt. blow) then
            Pbio(i) = (bag(i) - blow) / (bup - blow)
         else if (bag(i) .le. blow) then
            Pbio(i) = 0.0
         endif        
         Pbio(i) = max(0.0, min(Pbio(i),1.0))
250   continue

!-----------------------------------
! 2. Dependence on soil moisture, Pmoi

! Dryness factor for each soil layer

      do 300 j=1,nsoilay
         do 310 i=kpti,kptj
            wsoil(i,j) = wsoi(i,j) + wisoi(i,j)
            if (wsoil(i,j) .le. swilt(i,j)) then
               betad(i,j) = 0.0
            else if (wsoil(i,j) .gt. swilt(i,j) .and. wsoil(i,j) .lt. sfield(i,j)) then
               betad(i,j) = (wsoil(i,j) - swilt(i,j))
               betad(i,j) = betad(i,j) / (sfield(i,j) - swilt(i,j))
            else
               betad(i,j) = 1.0
            endif
            betad(i,j) = max(0.0, MIN(betad(i,j),1.0))
310      continue
300   continue

!   Soil moisture fire probability term, Pmoi
! Considers only 1st layer

      do 260 i=kpti,kptj
         Pmoi(i) = exp( -1.0 * PI * (betad(i,1) / betae)**2)
         Pmoi(i) = max(0.0, min(Pmoi(i),1.0))
260   continue

!----------------------------------------
! 3. dependence on ignition sources, Pign

! Ignition probability due to lightning
! Uses a random number

      do 270 i=kpti,kptj

! TODO replace Pligh by eq. (6), for now use random number
         Pligh = ran_uni( seedvec(i) )         

! Ignition probability by lightning + human

         Pign(i) = Pligh + (1 - Pligh) * Ph
         Pign(i) = max(0.0, min(Pign(i),1.0))

270   continue

!-----------------------
! 4. Total fire probability

      do 280 i=kpti,kptj
         Pfire(i) = Pbio(i) * Pmoi(i) * Pign(i)
         Ifire(i) = 1.0
280   continue

!-------------------------------------------------------
!   Estimation of burned area
!
! 1. Spread rate as funcion of wind speed, soil moisture
! 2. Length to breadth ratio of fire
! 3. Area burned in 1 day 
! 4. Average area burned (extinguishing prob)
!---------------------------------------------------------

!-----------------------
! 1. Spread rate

      do 430 i=kpti,kptj

         if (Ifire(i) .eq. 1.0) then

! Dependence of fire spread on root zone soil wetness

            if (betad(i,1) .gt. betae) then
               betas(i) = 1.0
            else
               betas(i) = betad(i,1) / betae
            endif

            hs(i) = (1.0 - betas(i))**2.0

!   Dependence of spread rate on wind speed
! Adjust value of alpha in params/fire

            wind(i) = ua(i) * 3.60     ! convert m/s to km/h
            ws(i) = 1.0 - ( (1.0 - g0) * exp(-1.0 * alpha * wind(i)**2) )

!   Spread rate of fire
! Adjust value of umax in params/fire

            srate(i) = umax * hs(i) * ws(i)

! Length to breadth ratio of fire

            lbratio(i) = 1.0 + 10.0 * (1.0 - exp(-0.017 * wind(i)))

!   Area burned in 1 day
! Considers t = 24 h

            abday(i) = (PI * 0.36 * 24 * 24 * srate(i)**2) / lbratio(i)

!   Total vegetated area
! Calculation of vegetated area accounts for upper and lower canopy
! fraction in a such way that don't favour any canopy

            sfrac(i) = fu(i) + fl(i)

            if (sfrac(i) .le. 1.0) then
               cfrac(i) = 1.0
            else
               cfrac(i) = 1.5 - (sfrac(i)/2)
            endif

            vegfrac(i) = cfrac(i) * sfrac(i)

            if (vegfrac(i) .gt. 1.0) then
               vegfrac(i) = 1.0
            endif

!   Average burned area
! depends on fire extinguishing probability
! Adjust exfire to an appropriate value in params/fire

            areax(i) = ((1.0 - exfire) * (2.0 - exfire)) / exfire**2.0

            ! TODO this still needs feed-back to vegetation and/or calibration/adjustments  
            ! because some points have ayburnfrac > 1
            burnfrac(i) = Pfire(i) * abday(i) * areax(i) * vegfrac(i) / reparea

         endif

!    Burned area can't be greater than vegetated area
! Verify if burned area is not greater than vegetated area

         if (burnfrac(i) .gt. vegfrac(i)) then
            burnfrac(i) = vegfrac(i)
         endif

430   continue

!   Burned area for each PFT
! Calculate for upper canopy (pfts 1 to 8) and lower (pfts 9 to 12)
! First lower canopy is burned
! If lower canopy area < burned area
! Upper canopy is burned

      do 520 i=kpti,kptj
         burnlow(i) = fl(i) * burnfrac(i)
         burnup(i) = burnfrac(i) - burnlow(i)
520   continue

      do 530 j=1,8
         do 531 i=kpti,kptj
            burnpft(i,j) = frac(i,j) * burnup(i)
531	 continue
530   continue

      do 540 j=9,12
         do 541 i=kpti,kptj
            burnpft(i,j) = frac(i,j) * burnlow(i)
541	 continue
540   continue

      do 550 j=1,12
         do 551 i=kpti,kptj
            sburnpft(i) = sburnpft(i) + burnpft(i,j)
551	 continue
550   continue

      do 560 i=kpti,kptj
         if (sburnpft(i) .gt. burnfrac(i)) then
            burnfrac(i) = sburnpft(i)
         endif
560   continue

      return

end subroutine firectem
