#include "inland_config.h"
!
! #    #  ######    ##     #####  #    #  ######  #####
! #    #  #        #  #      #    #    #  #       #    #
! #    #  #####   #    #     #    ######  #####   #    #
! # ## #  #       ######     #    #    #  #       #####
! ##  ##  #       #    #     #    #    #  #       #   #
! #    #  ######  #    #     #    #    #  ######  #    #
!
! TODO: We will probably need to support imthlast to denote last month for
!      restart, if iyrlast here denotes last year for restart as well.
! ---------------------------------------------------------------------
subroutine dailymon (seed, seed2,seed3,seed4, jdaily, iyrlast, nrun)
! ---------------------------------------------------------------------
!
! overview
!
! this routine generates daily weather conditions from monthly-mean
! climatic parameters
!
! specifically, this routine generates daily values of
!
!  - daily total precipitation
!  - daily maximum temperature
!  - daily minimum temperature
!  - daily average cloud cover
!  - daily average relative humidity
!  - daily average wind speed
!
! in order to generate daily weather conditions, the model uses a series
! of 'weather generator' approaches, which generate random combinations of
! weather conditions based upon the climatological conditions
!
! in general, this weather generator is based upon the so-called Richardson
! weather generator
!
! appropriate references include:
!
! Geng, S., F.W.T. Penning de Vries, and L. Supit, 1985:  A simple
! method for generating rainfall data, Agricultural and Forest
! Meteorology, 36, 363-376.
!
! Richardson, C. W. and Wright, D. A., 1984: WGEN: A model for
! generating daily weather variables: U. S. Department of
! Agriculture, Agricultural Research Service.
!
! Richardson, C., 1981: Stochastic simulation of daily
! precipitation, temperature, and solar radiation. Water Resources
! Research 17, 182-190.
!
! former common blocks, now they are fortran 90 modules

!
! FIX ME: Next version need merge subroutine dailymon with daily.
!

      use inland_parameters
      use inland_control, only: iyear, imonth, iday, iyear0
      use inland_combcs
      use inland_comatm
      use inland_comsum
      use inland_comveg
      use inland_comsatp
      use inland_subgrid
      use inland_comcrop
      use inland_comsoi
      use inland_com1d

      implicit none

! Arguments
      integer :: seed,     &
                 seed2,    &
                 seed3(32), &
                 seed4,     &
                 jdaily      ! 1 if reading in daily weather data
                             ! 0 if using random/statistical weather generator

! local variables
      integer :: it1w,     & ! indice of previous month (interpolation)
                 it2w,     & ! indice of following month (interpolation)
                 i,j,k,    & ! loop indice
                 ilpt        ! mlpt index

      real*8 :: rwork, aa, ab, tr1, tr2, s1, s2, s12, z, v, y, b3, b2, b1, x1, &
                amn,         & !
                omcloud,     & ! cloud cover
                omqd,        & ! humidity
                omtmax,      & ! maximum temperature
                dt,          & ! used for interpolation
                pwet,        & ! monthly-average probability of rainy day
                pwd,         & ! probability of a wet day after a dry day
                pww,         & ! probability of a wet day after a wet day
                rainpwd,     & ! average rainfall per wet day
                alpha,       & ! parameter for gamma function
                beta,        & ! parameter for gamma function
                tdm,         & !  mean daily mean temperature
                trngm,       & ! mean daily mean temperature
                tmaxm,       & ! mean maximum temperature
                tminm,       & ! mean minimum temperature
                tmaxd,       & ! maximum temperatures for dry days
                tmaxw,       & ! maximum temperatures for wet days
                tmaxe,       & ! 'expected' maximum temperature for today
                tmine,       & ! 'expected' minimum temperature for today
                tmaxs,       & ! standard deviation in maximum temperature (K)
                tmins,       & ! standard deviation in minimum temperature (K)
                cloudm,      & ! mean cloud cover for today (fraction)
                cloudd,      & ! dry day cloud cover
                cloudw,      & ! wet day cloud cover
                cloude,      & ! expected cloud cover today
                clouds,      & ! standard deviation of cloud fraction
                tdum,        & ! storage variable
                qdm,         & ! mean relative humidity
                qdd,         & ! dry day relative humidity
                qdw,         & ! wet day relative humidity
                qde,  & ! expected relative humidity (based on wet/dry decision)
                qdup,        & ! upper bound of humidity distribution function
                qdlow,       & ! lower bound of humidity distribution function
                eud        ! expected daily average wind speed from monthly mean

      real*8 :: a(3,3), b(3,3)
      real*8 :: ee(3), r(3), rr(3), x(3)
      integer :: iyrlast, nrun, iyrmon
      real*8 :: precipfac, dif
      real :: ran2,          & ! function random number generator
              rndnum,        & ! random number to decide if wet or dry day
              rn1,rn2,rn3,rn   ! random numbers

! define autocorrelation matrices for Richardson generator
!
! note that this matrix should be based upon a statistical
! analysis of regional weather patterns
!
! for global simulations, we use 'nominal' values
      data a /  0.600,  0.500,  0.005, &
                0.010,  0.250,  0.005, &
                0.020,  0.125,  0.250 /

      data b /  0.500,  0.250, -0.250, &
                0.000,  0.500,  0.250, &
                0.000,  0.000,  0.500 /

#define WGEN_COMSAT
#include "inland_comsat.h"

! ----------------------------------------------------------------------
! * * * initial setup for daily climate calculations * * *
! ----------------------------------------------------------------------
!
! define working variables
      rwork = (grav / rair / 0.0065)

! 'omega' parameters used to calculate differences in expected
! climatic parameters on wet and dry days
!
! following logic of weather generator used in the EPIC crop model
!
! omcloud -- cloud cover
! omqd    -- humidity
! omtmax  -- maximum temperature
      omcloud = 0.90    ! originally 0.90
      omqd    = 0.50    ! originally 0.50
      omtmax  = 0.75    ! originally 0.75

! calculate weighting factors used in interpolating climatological
! monthly-mean input values to daily-mean values
!
! this is a simple linear interpolation technique that takes into
! account the length of each month
      if (jdaily .eq. 0) then
         if (float(iday).lt.float(ndaypm(imonth)+1)/2.0) then
            it1w = imonth - 1
            it2w = imonth
            dt   = (float(iday) - 0.5) / ndaypm(imonth) + 0.5
         else
            it1w = imonth
            it2w = imonth + 1
            dt   = (float(iday) - 0.5) / ndaypm(imonth) - 0.5
         end if
         if (it1w.lt. 1) it1w = 12
         if (it2w.gt.12) it2w = 1
      end if

! initialize this year's values of gdd0, gdd5, tc, tw
      if (iday.eq.1 .and. imonth.eq.1) then
         ! we use 'mpt' in place of 'npoi' here: it's lend-lbeg+1!
         tcthis(:)=100.0
         twthis(:)=-100.0
         gdd0this(:)=0.0
         gdd5this(:)=0.0

         if(isimagro .gt. 0) then
            gdd0cthis(:)=0.0
            gdd8this(:)=0.0
            gdd10this(:)=0.0
            gdd11this(:)=0.0
            gdd12this(:)=0.0
!
! initialize variables to zero at beginning of year
! for crop types that do not grow over two calendar years
!
! constant for gdd based on 10 cm soil temperature (planting bed-
! used to calculate leaf emergence after planting
!
           do  i = lbeg, lend
              do  j = 1, npft
                 if (j .le. scpft-1) then  ! natural vegetation
                    ayanpp(i,j)     = 0.0
                 else if (j .ge. scpft) then
                    ayanpp(i,j)     = 0.0
                 endif
	          enddo
           enddo
         endif ! check for crop existence
      endif

     if(isimagro .gt. 0) then
!
! initialize this crop year's values
!
      do  i = lbeg, lend
         do  j = scpft, ecpft
            if (iday.eq.pcd(j).and.imonth.eq.pcm(j)) then
               if (exist(i,13).eq.1.and.j.eq.13) then
                  gddsoy(i,iyear-iyear0+5) = gddfzsoy(i)
                  consdays(i)=0
                  iniday(i)=9999
                  maxcons(i)=0
                  gsdays(i)=0
                  gddfzcorn(i)=0.0
                  gddfzsoy(i)=0.0
                  gddfzsgc(i)=0.0
                  gddfzplm(i)=0.0
               else if (exist(i,13).eq.0.and.j.eq.13) then
	          gddsoy(i,iyear-iyear0+5)=0.0
	       endif
!
	       if (exist(i,14).eq.1.and.j.eq.14) then
                  gddcorn(i,iyear-iyear0+5) = gddfzcorn(i)
                  consdays(i)=0
                  iniday(i)=9999
                  maxcons(i)=0
                  gsdays(i)=0
                  gddfzcorn(i)=0.0
                  gddfzsoy(i)=0.0
                  gddfzsgc(i)=0.0
                  gddfzplm(i)=0.0
               else if (exist(i,14).eq.0.and.j.eq.14) then
                  gddcorn(i,iyear-iyear0+5)=0.0
               endif
!
	       if (exist(i,16).eq.1.and.j.eq.16) then
                  gddsgc(i,iyear-iyear0+5)  = gddfzsgc(i)
                  consdays(i)=0
                  iniday(i)=9999
                  maxcons(i)=0
                  gsdays(i)=0
                  gddfzcorn(i)=0.0
                  gddfzsoy(i)=0.0
                  gddfzsgc(i)=0.0
                  gddfzplm(i)=0.0
               else if (exist(i,16).eq.0.and.j.eq.16) then
                  gddsgc(i,iyear-iyear0+5)=0.0
               endif

               if (exist(i,17).eq.1.and.j.eq.17) then
                        gddsgc(i,iyear-iyear0+5)  = gddfzsgc(i)
                        consdays(i)=0
                        iniday(i)=9999
                        maxcons(i)=0
                        gsdays(i)=0
                        gddfzcorn(i)=0.0
                        gddfzsoy(i)=0.0
                        gddfzsgc(i)=0.0
                        gddfzplm(i)=0.0
    else if (exist(i,17).eq.0.and.j.eq.17) then
                        gddplm(i,iyear-iyear0+5)=0.0
                     endif
!
               if (croplive(i,j) .eq. 0 ) then
                  gddplant(i,j)   = 0.0
                  gddtsoi(i,j)    = 0.0
               endif
            endif
         enddo
      enddo
!
    endif ! check for crop existence
! ----------------------------------------------------------------------
! * * * set daily climatic variables for entire domain * * *
! ----------------------------------------------------------------------

! if mlpt > 1 do for all gridcells
! else
! do for non-tile gridcells because of multiple calls to ran2 will generate diff. results
! and then replicate results to subgrid tiles
! in any case, loop from 1 to npoi1, because when mlpt = 1, npoi1 = npoi
! TODO fix this when lbeg != 1 and lend != npoi1
      do 200 i = 1, npoi1

! ----------------------------------------------------------------------
! * * * use weather generator to create daily statistics * * *
! ----------------------------------------------------------------------
         if (jdaily .eq. 0) then

! ----------------------------------------------------------------------
! (1) determine if today will rain or not (following Geng et al.)
! ----------------------------------------------------------------------
!
! implement simple first-order Markov-chain precipitation generator logic
! based on Geng et al. (1986), Richardson and Wright (1984),
! and Richardson (1981)
!
! basically, this allows for the probability of today being a wet day
! (a day with measureable precipitation) to be a function of what
! yesterday was (wet or dry)
!
! the logic here is that it is more likely that a wet day will follow
! another wet day -- allowing for 'storm events' to persist
!
! calculate monthly-average probability of rainy day
            pwet = max (dble(1.), xinwetmon(i,imonth)) / ndaypm(imonth)

! estimate the probability of a wet day after a dry day
            pwd = 0.75 * pwet

! estimate the probability of a wet day after a wet day
            pww = 0.25 + pwd

! Beginning of block of code that figures out daily precip for
! entire month on the first day of each month
            if (iday .eq. 1) then

! Verify the dataset consistency especially when using interannual anomalies
! of precipitations (negative values, or too few wet days in a rainy month)
               xinprecmon(i, imonth) = max(dble(0.), xinprecmon(i, imonth))
               xinwetmon(i, imonth) = max(dble(1.), xinwetmon(i, imonth))

9000           continue

! Initialize monthly parameters back to zero
               iwetdaysum(i) = 0
               precipdaysum(i) = 0

               do 210 j = 1, 31
                  iwetday(i,j) = 0
                  precipday(i,j) = 0
210            continue

! Loop through number of days in this month and determine precip
               do 220 j = 1, ndaypm(imonth)

! decide if today is a wet day or a dry day using a random number
                  rndnum = ran2(seed,seed2,seed3,seed4)

! If it is the first day of the month do not look at previous day
                  if (j .eq. 1) then
                     if (dble(rndnum) .le. pwd) then
                        iwetday(i,j) = 1
                        iwetdaysum(i) = iwetdaysum(i) + 1
                     else
                        iwetday(i,j) = 0
                     endif
                  else

! If it is not the first day, look at yesterday's wet/dry index to help
! determine if today is wet or dry
                     if (iwetday(i,j-1) .eq. 0) then
                        if (dble(rndnum).le.pwd) then
                           iwetday(i,j) = 1
                           iwetdaysum(i) = iwetdaysum(i) + 1
                        endif
                     else
                        if (dble(rndnum).gt.pww) iwetday(i,j) = 0
                     endif
                  endif

! ----------------------------------------------------------------------
! (2) determine today's precipitation amount (following Geng et al.)
! ----------------------------------------------------------------------
!
! if it is going to rain today
                  if (iwetday(i,j) .eq. 1) then

! calculate average rainfall per wet day
                     rainpwd = xinprecmon(i,imonth) * ndaypm(imonth) / &
                               max (dble(0.1), xinwetmon(i,imonth))

! randomly select a daily rainfall amount from a probability density
! function of rainfall
!
! method i --
!
! use the following technique from Geng et al. and Richardson
! to distribute rainfall probabilities
!
! pick a random rainfall amount from a two-parameter gamma function
! distribution function
!
! estimate two parameters for gamma function (following Geng et al.)
                     beta  = max (dble(1.0), -2.16 + 1.83 * rainpwd)
                     alpha = rainpwd / beta

! determine daily precipitation amount from gamma distribution function
! (following WGEN code of Richardson and Wright (1984))
                     aa = 1.0 / alpha
                     ab = 1.0 / (1.0 - alpha)
                     tr1 = exp(-18.42 / aa)
                     tr2 = exp(-18.42 / ab)

12                   rn1 = ran2(seed,seed2,seed3,seed4)
                     rn2 = ran2(seed,seed2,seed3,seed4)

! CD: rewrote parts of prehistoric code in fortran 77
                     if ((dble(rn1) - tr1) .le. 0) then
                        s1 = 0.0
                     else
                        s1 = dble(rn1)**aa
                     end if

                     if ((dble(rn2) - tr2) .le. 0) then
                        s2 = 0.0
                     else
                        s2 = dble(rn2)**ab
                     end if

!                 if (rn1 - tr1) 61, 61, 62
! 61              s1 = 0.0
!                 go to 63
! 62              s1 = rn1**aa
! 63              if (rn2 - tr2) 64, 64, 65
! 64              s2 = 0.0
!                 go to 66
! 65              s2 = rn2**ab
!
! 66                   s12 = s1 + s2
                     s12 = s1 + s2
!
                     if (s12 - 1.0)  13, 13, 12
13                      z = s1 / s12
                        rn3 = ran2(seed,seed2,seed3,seed4)
                        precipday(i,j) = -z * dble(log(rn3)) * beta

! method ii --
!
! here we use a one-parameter Weibull distribution function
! following the analysis of Selker and Haith (1990)
!
! Selker, J.S. and D.A. Haith, 1990: Development and testing of single-
! parameter precipitation distributions, Water Resources Research,
! 11, 2733-2740.
!
! this technique seems to have a significant advantage over other
! means of generating rainfall distribution functions
!
! by calibrating the Weibull function to U.S. precipitation records,
! Selker and Haith were able to establish the following relationship
!
! the cumulative probability of rainfall intensity x is given as:
!
! f(x) = 1.0 - exp(-(1.191 x / rainpwd)**0.75)
!
! where x       : rainfall intensity
!       rainpwd : rainfall per wet day
!
! using transformation method, take uniform deviate and convert it to a
! random number weighted by the following Weibull function
!
!          rndnum = ran2(seed,seed2,seed3,seed4)
!
!          precip(i) = rainpwd / 1.191 * (-log(1.0 - rndnum))**1.333333
!
! bound daily precipitation to "realistic" range
!
! lower end is determined by definition of a 'wet day' (at least
! 0.25 mm of total precipitation)
!
! upper end is to prevent inland from blowing up
                     precipday(i,j) = max (precipday(i,j),dble(0.25))    ! min =   0.25 mm/day
                     precipday(i,j) = min (precipday(i,j),dble(150.00))  ! max = 150.00 mm/day
!
! Back to beginning of month loop, this is the end of it
                  endif

! Add today's precip to the monthly summation
                  precipdaysum(i) = precipdaysum(i) + precipday(i,j)

220            continue
!
! Adjust daily precip amounts (using precipfac) so that the monthly
! summation equals the input precip amount, when using interannual
! anomalies
               if ((precipdaysum(i) .eq. 0).AND.(xinprecmon(i,imonth) .gt. 0)) then
                  rndnum = 1.0 + (float(ndaypm(imonth)) - 1.0) * ran2(seed,seed2,seed3,seed4)
                  iwetday(i,nint(rndnum)) = 1
                  precipday(i,nint(rndnum)) = xinprecmon(i,imonth) *       &
                                              float(ndaypm(imonth))
                  precipdaysum(i) = precipday(i,nint(rndnum))
                  iwetdaysum(i) = 1
               end if

               precipfac = (xinprecmon(i,imonth)*float(ndaypm(imonth)))  / &
                           max(dble(0.01),precipdaysum(i))

!   The test below determines if there were more than 365 rainy days for the
! year. If so, it goes back and reestimate a better number.
!   In the case the input data map has a NaN on a grid point that should have
! a valid value, it will incur in an infinite loop as there will be a negative
! number of rainy days, or a too great number that does not change every '9000
! loop'. Thus, now we check if such a situation happens and stop the model
! warning the user of the fact. More details like the coordinate of  the point
! could help identify the corrupt data on the precipitation input file, but we
! are keeping just safe for now. The affected 'npoi' is printed. From this
! number the user can infer the gridpoint global coordinates and locate it on
! the data files should this problem ever happens again.
! - fabricio 20110826
               do 230 j=1,ndaypm(imonth)
                  precipday(i,j) = precipday(i,j) * precipfac
                  if (precipday(i,j).gt.360) then
                     if (xinwetmon(i,imonth) .lt. ndaypm(imonth)) then
                        xinwetmon(i,imonth) = xinwetmon(i,imonth) + 1
                        pwet = xinwet(i,imonth) / ndaypm(imonth)
                        pwd = 0.75 * pwet
                        pww = 0.25 + pwd
                        print *, 'WARNING: goto 9000a', i,             &
                                 int(xinwetmon(i,imonth)), iwetdaysum(i), &
                                 int(precipday(i,j))
                          if(precipday(i,j).le.0. .or. &
                             precipday(i,j).ge.100000.) then
                             print *,"xinprecmon(i,month):",xinprecmon(i,imonth)," xinwetmon():", xinwetmon(i,imonth)," i:", i, " month:",imonth
                             print *,"Unexpected value from precipitation input data."
                             stop 1
                          endif
                        goto 9000
                     else
                        print *, 'WARNING: goto 9000b', i,             &
                                 int(xinwetmon(i,imonth)), iwetdaysum(i), &
                                 int(precipday(i,j))
                          if(precipday(i,j).le.0. .or. &
                             precipday(i,j).ge.100000.) then
                             print *,"xinprecmon(i,month):",xinprecmon(i,imonth)," xinwetmon():", xinwetmon(i,imonth)," i:", i, " month:",imonth
                             print *,"Unexpected value from precipitation input data."
                             stop 2
                          endif
                        goto 9000
                     end if
                  end if
230            continue

! Verification of the weather generator algorithm
               iwetdaysum(i) = 0
               precipdaysum(i) = 0.

               do 240 j=1,ndaypm(imonth)
                  precipdaysum(i) = precipdaysum(i) + precipday(i,j)
                  iwetdaysum(i) = iwetdaysum(i) + iwetday(i,j)
240            continue

               dif = precipdaysum(i) - xinprecmon(i,imonth) * float(ndaypm(imonth))

               if ((dif.lt.-0.1).or.(dif.gt.0.1)) then
                  print *, 'ERROR in DAILY:', i, precipdaysum(i),    &
                           xinprecmon(i,imonth)* float(ndaypm(imonth)), &
                           iwetdaysum(i), xinwetmon(i,imonth)
               end if
! end of the verification
            end if               !end of the iday loop

! Relate today's iwetday and precipday to iwet and precip that will be
! used below
            iwet(i) = iwetday(i,iday)
            precip(i) = precipday(i,iday)

! ----------------------------------------------------------------------
! (3) estimate expected minimum and maximum temperatures
! ----------------------------------------------------------------------
!
! first determine the expected maximum and minimum temperatures
! (from climatological means) for this day of the year
!
! mean daily mean temperature (K)

            if (iyear.eq.iyrmon) then
                 if(imonth.eq.1) then
                 xintmon(i,it1w) = xint(i,it1w)
                 xintrngmon(i,it1w) = xintrng(i,it1w)
                 endif
            endif
            if (iyear.eq.iyrlast + nrun) then
                if(imonth.eq.12) then
                 xintmon(i,it2w) =  xint(i,it2w)
                 xintrngmon(i,it2w) = xintrng(i,it2w)
                endif
            endif

            tdm = xintmon(i,it1w) + dt * (xintmon(i,it2w) - xintmon(i,it1w)) + 273.16

! mean daily temperature range (K)
            trngm = xintrngmon(i,it1w) + dt * (xintrngmon(i,it2w) - xintrngmon(i,it1w))

! mean minimum and maximum temperatures
            tmaxm = tdm + 0.56 * trngm
            tminm = tdm - 0.44 * trngm

! modify maximum temperatures for wet and dry days
            if (pwet .ne. 0.0) then
               tmaxd = tmaxm + pwet * omtmax * trngm
               tmaxw = tmaxd -        omtmax * trngm
            else
               tmaxd = tmaxm
               tmaxw = tmaxm
            endif

! set the 'expected' maximum and minimum temperatures for today
!
! note that the expected minimum temperatures are the same for
! both wet and dry days
            if (iwet(i).eq.0) tmaxe = tmaxd
            if (iwet(i).eq.1) tmaxe = tmaxw

            tmine = tminm

! estimate variability in minimum and maximum temperatures
!
! tmaxs : standard deviation in maximum temperature (K)
! tmins : standard deviation in minimum temperature (K)
!
! Regression is based on analysis of 2-m air temperature data from the
! NCEP/NCAR reanalysis (1958-1997) for 294 land points over central
! North America (24N-52N, 130W-60W, 0.5-degree resolution): Daily max
! and min temperatures were calculated for each land point from daily
! mean temperature and temperature range. Anomalies were calculated
! by subtracting similar max and min temperatures calculated from
! monthly mean temperature and range (interpolated to daily). The 40
! years of anomalies were then binned by month and the standard
! deviation calculated for each month. The 294 x 12 standard
! deviations were then regressed against the 3528 long-term monthly
! mean temperatures.
!
! note: the values are bound to be greater than 1.0 K
! (at the very least they must be bound so they don't go below zero)
            tmaxs = max (dble(1.0), -0.0713 * (tdm - 273.16) + 4.89)
            tmins = max (dble(1.0), -0.1280 * (tdm - 273.16) + 5.73)

! ----------------------------------------------------------------------
! (4) estimate expected cloud cover
! ----------------------------------------------------------------------
!
! the formulation of dry and wet cloud cover has been
! derived from the weather generator used in the epic crop model
!
! cloudm : mean cloud cover for today
! cloudd : dry day cloud cover
! cloudw : wet day cloud cover
! cloude : expected cloud cover today
!
! Verify the data set consistency when using interannual anomalies of
! cloudiness (values under 0 % or over 100 %)
            if (iday.eq.1) then
               xincld(i,it1w) = max (dble(0.), xincld(i,it1w))
               xincld(i,it1w) = min (dble(100.), xincld(i,it1w))
               xincld(i,it2w) = max (dble(0.), xincld(i,it2w))
               xincld(i,it2w) = min (dble(100.), xincld(i,it2w))
               xincldmon(i,imonth) = max (dble(0.), xincldmon(i,imonth))
               xincldmon(i,imonth) = min (dble(100.), xincldmon(i,imonth))
            end if

            if (iyear.eq.iyrmon) then
                 if(imonth.eq.1) then
                 xincldmon(i,it1w) = xincld(i,it1w)
                 endif
            endif
            if (iyear.eq.iyrlast + nrun) then
                if(imonth.eq.12) then
                 xincldmon(i,it2w) =  xincld(i,it2w)
                endif
            endif

! monthly mean cloud cover (%)
            cloudm = xincldmon(i,it1w) + dt * (xincldmon(i,it2w) - xincldmon(i,it1w))

! convert from percent to fraction
            cloudm = cloudm / 100.0

! adjust cloud cover depending on dry day / rainy day
! following logic of the EPIC weather generator code
            if (pwet .ne. 0.0) then
               cloudd = (cloudm - pwet * omcloud) / (1.0 - pwet * omcloud)
               cloudd = min (dble(1.0), max (dble(0.0), cloudd))
               cloudw = (cloudm - (1.0 - pwet) * cloudd) / pwet
            else
               cloudd = cloudm
               cloudw = cloudm
            endif
            if (iwet(i).eq.0) cloude = cloudd
            if (iwet(i).eq.1) cloude = cloudw

! estimate variability in cloud cover for wet and dry days
! following numbers proposed by Richardson
!
! clouds : standard deviation of cloud fraction
            if (iwet(i).eq.0) clouds = 0.24 * cloude
            if (iwet(i).eq.1) clouds = 0.48 * cloude

! ----------------------------------------------------------------------
! (5) determine today's temperatures and cloud cover using
!     first-order serial autoregressive technique
! ----------------------------------------------------------------------
!
! use the Richardson (1981) weather generator approach to simulate the
! daily values of minimum / maximum temperature and cloud cover
!
! following the implementation of the Richardson WGEN weather generator
! used in the EPIC crop model
!
! this approach uses a multivariate generator, which assumes that the
! perturbation of minimum / maximum temperature and cloud cover are
! normally distributed and that the serial correlation of each
! variable may be described by a first-order autoregressive model
!
! generate standard deviates for weather generator
            do 111 j = 1, 3
 31            rn1 = ran2(seed,seed2,seed3,seed4)
               rn2 = ran2(seed,seed2,seed3,seed4)
               v = sqrt (-2.0 * dble(log(rn1))) * dble(cos(6.283185 * rn2))
               if (abs(v) .gt. 2.5) go to 31
               ee(j) = v
111         continue

! zero out vectors
            do 121 j = 1, 3
               r(j)  = 0.0
               rr(j) = 0.0
121         continue

! update working vectors
            do 131 j = 1, 3
               do 141 k = 1, 3
                  r(j)  = r(j)  + b(j,k) * ee(j)
                  rr(j) = rr(j) + a(j,k) * xstore(i,k)
141            continue
131         continue
!
! solve for x() perturbation vector and save current vector
! into the xim1() storage vector (saved for each point)
            do 151 j = 1, 3
               x(j) = r(j) + rr(j)
               xstore(i,j) = x(j)
151         continue

! determine today's minimum and maximum temperature
            tmax(i)  = tmaxe + tmaxs * x(1)
            tmin(i)  = tmine + tmins * x(2)

! if tmin > tmax, then switch the two around
            if (tmin(i).gt.tmax(i)) then
               tdum    = tmax(i)
               tmax(i) = tmin(i)
               tmin(i) = tdum
            endif

! daily average temperature
            td(i) = 0.44 * tmax(i) + 0.56 * tmin(i)

! determine today's cloud cover
            cloud(i) = cloude + clouds * x(3)

! constrain cloud cover to be between 0 and 100%
            cloud(i) = max (dble(0.0), min (dble(1.0), cloud(i)))

! ----------------------------------------------------------------------
! (6) estimate today's surface atmospheric pressure
! ----------------------------------------------------------------------
!
! simply a function of the daily average temperature and topographic
! height -- nothing fancy here
            psurf(i) = 101325.0 * ((td(i) - 0.0065 * xintopo(i))/td(i))**rwork

! ----------------------------------------------------------------------
! (7) estimate today's relative humidity
! ----------------------------------------------------------------------
!
! the formulation of dry and wet relative humidities has been
! derived from the weather generator used in the epic crop model
!
! qdm : mean relative humidity
! qdd : dry day relative humidity
! qdw : rainy day relative humidity
! qde : expected relative humidity (based on wet/dry decision)
!
! Verify the data set consistency when using interannual anomalies of
! relative humidity (values over 100 % or under 0 %)
            if (iday.eq.1) then
               xinq(i,it1w) = max (dble(0.), xinq(i,it1w))
               xinq(i,it1w) = min (dble(100.), xinq(i,it1w))
               xinq(i,it2w) = max (dble(0.), xinq(i,it2w))
               xinq(i,it2w) = min (dble(100.), xinq(i,it2w))
               xinqmon(i,imonth) = max (dble(0.), xinqmon(i,imonth))
               xinqmon(i,imonth) = min (dble(100.), xinqmon(i,imonth))
            end if

            if (iyear.eq.iyrmon) then
                 if(imonth.eq.1) then
                    xinqmon(i,it1w) = xinq(i,it1w)
                 endif
            endif

            if (iyear.eq.iyrlast + nrun) then
                if(imonth.eq.12) then
                   xinqmon(i,it2w) =  xinq(i,it2w)
                endif
            endif

! mean relative humidity (%)
            qdm = xinqmon(i,it1w) + dt * (xinqmon(i,it2w) - xinqmon(i,it1w))

! convert from percent to fraction
            qdm = qdm / 100.0

! adjust humidity depending on dry day / rainy day
! following logic of the EPIC weather generator code
            if (pwet .ne. 0.0) then
               qdd = (qdm - pwet * omqd) / (1.0 - pwet * omqd)
               if (qdd .lt. 0.2) then
                  qdd = 0.2
                  if (qdd .gt. qdm) qdm = qdd
               endif
               qdd = min(dble(1.0), qdd)
               qdw = (qdm - (1.0 - pwet) * qdd) / pwet
            else
               qdd = qdm
               qdw = qdm
            endif

            if (iwet(i).eq.0) qde = qdd
            if (iwet(i).eq.1) qde = qdw

! estimate lower and upper bounds of humidity distribution function
! following logic of the EPIC weather generator code
            qdup  = qde + (1.0 - qde) * exp (qde - 1.0)
            qdlow = qde * (1.0 - exp (-qde))

! randomly select humidity from triangular distribution function
! following logic of the EPIC weather generator code
            rn = ran2(seed,seed2,seed3,seed4)
            y  = 2.0 / (qdup - qdlow)
            b3 = qde  - qdlow
            b2 = qdup - qde
            b1 = dble(rn) / y
            x1 = y * b3 / 2.0

            if (dble(rn).gt.x1) then
               qd(i) = qdup  - sqrt (b2 * b2 - 2.0 * b2 * (b1 - 0.5 * b3))
            else
               qd(i) = qdlow + sqrt (2.0 * b1 * b3)
            endif

! adjust daily humidity to conserve monthly mean values
!
! note that this adjustment sometimes gives rise to humidity
! values greater than 1.0 -- which is corrected below
            amn = (qdup + qde + qdlow) / 3.0
            qd(i) = qd(i) * qde / amn

! constrain daily average relative humidity
            qd(i) = max (dble(0.30), qd(i))
            qd(i) = min (dble(0.99), qd(i))

! convert from relative humidity to specific humidity at
! daily mean temperature
            qd(i) = qd(i) * qsat(esat(td(i)), psurf(i))

! ----------------------------------------------------------------------
! (8) estimate today's daily average wind speed
! ----------------------------------------------------------------------
!
! first estimate the expected daily average wind speed (from monthly
! means)

            if (iyear.eq.iyrmon) then
                 if(imonth.eq.1) then
                 xinwindmon(i,it1w) = xinwind(i,it1w)
                 endif
            endif
            if (iyear.eq.iyrlast + nrun) then
                if(imonth.eq.12) then
                 xinwindmon(i,it2w) =  xinwind(i,it2w)
                endif
            endif

            eud = xinwindmon(i,it1w) + dt * (xinwindmon(i,it2w) - xinwindmon(i,it1w))

! following logic of the EPIC weather generator
! select random wind speed following this equation
            ud(i) = 1.13989 * eud * (-log(ran2(seed,seed2,seed3,seed4)))**0.30

! constrain daily wind speeds to be between 2.5 and 10.0 m/sec
            ud(i) = max (dble(2.5), min (dble(10.0), ud(i)))

! ----------------------------------------------------------------------
! * * * use real daily climate data * * *
! ----------------------------------------------------------------------
         else

! use basic daily climate data, converting units
!
! daily total precipitation
!
! Here we multiply xinprecd, the daily fraction of precip calculated from
! the NCEP dataset, by xinprec, the total monthly amount of precip taken from
! the CRU05 dataset to obtain our derived daily precip amount. Also do a check
! to see if the daily precip exceeds 360mm (as is done in the daily weather
! generator) ... no correction is made, only a warning is printed
!            precip(i) = (xinprec(i,imonth) * ndaypm(imonth)) * xinprecd(i)
            precip(i) = xinprecd(i)
!           if (precip(i) .gt. 360) then
!              print *, 'WARNING: daily precip exceeds 360mm for'
!              print *, 'year, month, day, gridcell = '
!              print *, iyear, imonth, iday, i
!           endif
!
! daily average temperatures
!
! Here we add the NCEP temperature anomaly to the CRU05 monthly anomaly
! The trange NCEP anomaly must also be multiplied by the climatological
! CRU05 trange in order to determine tmax and tmin
            !td(i) = xint(i,imonth) + 273.16 + xintd(i)
            td(i) = xintd(i) + 273.16
	    !trngm = min (dble(44.0), (xintrng(i,imonth) * xintrngd(i)))
	    !tmax(i) = td(i) + 0.56 * trngm
            !tmin(i) = td(i) - 0.44 * trngm
            tmax(i) = xintmaxd(i) + 273.16
            tmin(i) = xintmind(i) + 273.16

! daily average cloud cover
!
! Here we add the NCEP cloud anomaly to the monthly anomaly from CRU05
! before converting percentage of cover to fraction
! We also bound cloud cover fraction between 0 and 1
            !cloud(i) = (xincld(i,imonth) + xincldd(i)) * 0.01
            cloud(i) = xincldd(i) *0.01
            cloud(i) = min (cloud(i), dble(1.0))
            cloud(i) = max (dble(0.0), cloud(i))

! compute surface atmospheric pressure
            psurf(i) = 101325.0 * ((td(i) - 0.0065 * xintopo(i))/td(i))**rwork


! daily average specific humidity
!
! First we must convert relative humidity to a fraction and then convert
! the fraction to specific humidity
! Then we can multiply the NCEP daily anomaly by the CRU05 monthly anomaly
            !humidfrac = xinq(i,imonth) / 100.
            qd(i) = xinqd(i) /100.
            qd(i) = qd(i) * qsat(esat(td(i)),psurf(i))

! daily average wind speed
!
! Here we multiply the NCEP fraction of windspeed by the CRU05
! climatological monthly windspeed
            !ud(i) = xinwind(i,imonth) * xinwindd(i)
            ud(i) = xinwindd(i)
            ud(i) = max (dble(2.5), min (dble(10.0), ud(i)))

!

! TET change 10/29/02
! Input windspeed is assumed to be measured at canopy height
! (ztop(i,2))
! IBIS now has siga = 0.991 which is ~93m height.
! Adjust windspeed according to log-wind profile.
! So let
! displacement height, d = 0.65(ztop(i,2))
! roughness length, zm = 0.1(ztop(i,2))
! Equation found in Bonan (2002) and other texts:
! u(z2)-u(z1)=(frictionvelocity/k)*log((z2-d)/(z1-d))
! Use log-wind to solve for frictionvelocity and substitute
! into above equation:
! u(z2) = u(z1)*(1.+(log((z2-d)/(z1-d))) / (log((z1-d)/zm)))
! Use canopyheight = z1 = ztop(i,2), and z2=93m
! and substitute d and zm to get equation dependent on canopy height:

        if(isimagro .gt. 0) then
           ud(i) = ud(i)*(1. + 0.79824*log((za(i)-0.65*ztop(i,2))/ &
                        (0.35*ztop(i,2))))
        endif

!
! Decided to use a canopy height of 20m over US after finding
! references of mature forest height between 10m and 35m (we are
! trying to simulate forests before human interference)
! with ztop(i,2) = 20m, the input windspeed is adjusted by a
! factor of 3.21
!
         end if

     if(isimagro .gt. 0)then

! * * *
      if (iyear.ge.istyear .and. iyear.le.istend ) then

      if(stinprecd(i).ge.0.) precip(i) = stinprecd(i)

      if(stintmax(i).ge.-40) tmax(i) = stintmax(i)
      if(stintmin(i).ge.-40) tmin(i) = stintmax(i)
      if(stintd(i).ge.-40) td(i) = stintd(i)


      if(stincldd(i).ge.0 .and.stincldd(i).le.1) &
        cloud(i) = stincldd(i)

      if(stintd(i).ge.-40) &
           psurf(i) = 101325.0 * &
                  (stintd(i) / (td(i) + 0.0065 * xintopo(i))) ** rwork

      if(stinqd(i).ge.0.and.stintd(i).ge.-40) then

      qdm = stinqd(i)


! convert from percent to fraction
          qdm = qdm / 100.0

! constrain daily average relative humidity
          qd(i) = max (0.30, qd(i))
          qd(i) = min (0.99, qd(i))
	print*,'warning'

!
! convert from relative humidity to specific humidity at
! daily mean temperature
          qd(i) = qd(i) * qsat(esat(stintd(i)), psurf(i))
!
	endif
       if(stinwindd(i).ge.0.) then

         ud(i) = stinwindd(i)

          ud(i) = ud(i)*(1. + 0.79824*log((93.-0.65*ztop(i,2))/ &
                  (0.35*ztop(i,2))))

	endif ! for wind

	stinprecd(i) = -999.0
	stintd(i)    = -999.0
	stintmax(i)  = -999.0
	stintmin(i)  = -999.0
	stincldd(i)  = -999.0
	stinqd(i)    = -999.0
	stinwindd(i) = -999.0
      endif
    endif ! check for crop existence
!
! ----------------------------------------------------------------------
! * * * other daily climate calculations * * *
! ----------------------------------------------------------------------
!
! calculated temperature extremes -- for vegetation limits (deg c)
!
! for this purpose, use the 10-day running mean temperature
         tcthis(i) = min (tcthis(i), (a10td(i) - 273.16))
         twthis(i) = max (twthis(i), (a10td(i) - 273.16))

! update this year's growing degree days
         gdd0this(i) = gdd0this(i) + max(dble(0.), (td(i) - 273.16))
         gdd5this(i) = gdd5this(i) + max(dble(0.), (td(i) - 278.16))

       if(isimagro .gt. 0) then

        gdd0cthis(i) = gdd0cthis(i) + max (0.0 ,(td(i) - baset(15)))    ! wheat
        gdd8this(i)  = gdd8this(i)  + max (0.0 ,(td(i) - baset(14)))    ! maize
        gdd10this(i) = gdd10this(i) + max (0.0 ,(td(i) - baset(13)))    ! soybean
        gdd11this(i) = gdd11this(i) + max (0.0 ,(td(i) - baset(17)))    ! oil palm
        gdd12this(i) = gdd12this(i) + max (0.0 ,(td(i) - baset(16)))    ! sugarcane

!
! form calculations of number of growing degree days between frost events
!
! events (e.g., tmin(i) .le. -2.2 C) this is applicable to CRU data
! differences exist when using a combination of CRU/NCEP
! -2.2 used as a threshold from Schwartz and Reiter, 2000.  International
! Journal of Climatology, 20: 929-932.  Changes in North American Spring
!
!
        if (tmin(i) .ge. 273.16) then
          consdays(i) = consdays(i) + 1
          maxcons(i)  = max(maxcons(i), consdays(i))
          if (maxcons(i) .eq. consdays(i)) then
!
            iniday(i) = cdays(i)+1 - maxcons(i)
          endif
!
          daygddc(1,cdays(i)) = min(30.0,max(0.0, (td(i) - baset(14))))
          daygdds(1,cdays(i)) = min(30.0,max(0.0, (td(i) - baset(13))))
          daygddsgc(1,cdays(i)) = min(30.0,max(0.0, (td(i) - baset(16))))

!
         else
            consdays(i) = 0
         endif
!
         if (cdays(i).eq.365) then
!Growing-degree day
            if (iniday(i) .eq. 9999) then
               iniday(i) = cdays (i)
               maxcons(i) = 1
            elseif (iniday(i) .eq. 0) then
	       iniday(i)=1
            endif
               endday(i) = iniday(i) + maxcons(i)-1

         do 125 k = iniday(i), endday(i)

             gddfzcorn(i) =  gddfzcorn(i) + daygddc(1,cdays(i))
             gddfzsoy(i)  =  gddfzsoy(i)  + daygdds(1,cdays(i))
             gddfzsgc(i) =  gddfzsgc(i)   + daygddsgc(1,cdays(i))
             gsdays(i) = gsdays(i) + 1


 125      continue

              elseif (cdays(i).eq.366) then

             if (iniday(i).eq.0) then
                 iniday(i) = 1
             endif

             endday(i) = iniday(i) + maxcons(i)-1

             gddfzcorn(i) = 0.
             gddfzsoy(i)  = 0.
             gddfzsgc(i) =  0.
             gsdays(i) =    0.



         do  k = iniday(i), endday(i)

             gddfzcorn(i) =  gddfzcorn(i) + daygddc(1,cdays(i))
             gddfzsoy(i)  =  gddfzsoy(i)  + daygdds(1,cdays(i))
             gddfzsgc(i) =  gddfzsgc(i)   + daygddsgc(1,cdays(i))
             gsdays(i) = gsdays(i) + 1
!
!
          enddo
!
        endif
! accumulate growing degree days for planted crops past planting
!
        do 150 j = scpft, ecpft
!
! for crops except winter wheat
! be careful with rotations
!

       if (croplive(i,j).eq.1.0.and.((j.eq.13.or.j.eq. 14.or.j.eq. 16) &
             .or.    (iwheattype .eq. 1 .and. j .eq. 15))) then
!
             gddplant(i,j) = gddplant(i,j) + max(0.0, min(td(i) &
                                           - baset(j), mxtmp(j)))
!
             gddtsoi(i,j)  = gddtsoi(i,j) + max(0.0, min(tsoi(i,1) &
                                          - baset(j), mxtmp(j)))

!
! if winter wheat is planted, reduce thermal time accumulation
! by vernalization factor (calculated in crops.f) until crop
! is fully vernalized (cold-hardened)
!
          else if (croplive(i,j) .eq. 1.0 .and. j .eq. 15 .and. &
                   iwheattype .eq. 2) then
!
             gddplant(i,j) = gddplant(i,j) + vf(i) * max(0.0, min(td(i) &
                                           - baset(j), mxtmp(j)))
!
             gddtsoi(i,j)  = gddtsoi(i,j)  +  vf(i) * max(0.0, min(tsoi(i,1) &
                                           - baset(j), mxtmp(j)))
!
        endif

	if(j.eq.16) then

	  ik(i)=ik(i)+1
	 if(cropy(i).eq.1.and.idpp(i,j).eq.1) then
!	 print*,'GDD_start',i,iyear,jday
         gddpl15(i)=0
	 ik(i)=1
	endif

       if (ik(i).le.(mxmat(16)-30)) then
       gddpl15(i)= gddpl15(i) + max(0.0, min(td(i)- baset(j), mxtmp(j)))
	endif

       if (ik(i).eq.(mxmat(16)-30)) then
!	print*,'GDD_1/5 anos',i,iyear,jday,gddsgcp(i,1),gddpl15(i),mxmat(16)-30
       gddsgcp(i,1)= (gddsgcp(i,1)+ gddpl15(i))/2
!	print*,'GDD_1/5 anos',i,iyear,jday,gddsgcp(i,1)
	endif

	endif

!
 150  continue
!c
!sant - out for geovane
!	if(i.eq.1)write(42,*)iyear,jday,xinprec(i,imonth),precip(i), &
!      xint(i,imonth), td(i)-273.18

   endif ! check for crop existence

 200  continue

! if using subgrid tiles, replicate results to subgrid tiles
if ( mlpt .gt. 1 ) then
! TODO fix this when lbeg != 1 and lend != npoi1
      do i = 1, npoi1
         do ilpt = 2,mlpt
            j = subgrid_get_index(i,ilpt)
            if ( j .ne. 0 ) then

               ! not sure these vars are needed
               iwetdaysum(j) = iwetdaysum(i)
               precipdaysum(j) = precipdaysum(i)
               iwetday(j,:) = iwetday(i,:)
               precipday(j,:) = precipday(i,:)
               iwet(j) = iwet(i)

               precip(j) = precip(i)
               tmax(j) = tmax(i)
               tmin(j) = tmin(i)
               td(j) = td(i)
               cloud(j) = cloud(i)
               psurf(j) = psurf(i)
               qd(j) = qd(i)
               ud(j) = ud(i)
               tcthis(j) = tcthis(i)
               twthis(j) = twthis(i)
               gdd0this(j) = gdd0this(i)
               gdd5this(j) = gdd5this(i)
               gdd8this(j) = gdd8this(i)
               gdd10this(j) = gdd10this(i)
               gdd11this(j) = gdd11this(i)
               gdd12this(j) = gdd12this(i)
               gddplant(j,:) = gddplant(i,:)
               gddtsoi(j,:) = gddtsoi(i,:)
               gddpl15(j) = gddpl15(i)
               gddsgcp(j,:) = gddsgcp(i,:)
            end if
         end do
      end do
endif
! return to main program
      return
end subroutine dailymon
