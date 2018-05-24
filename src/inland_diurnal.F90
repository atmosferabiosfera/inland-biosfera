#include "inland_config.h"
!
! #    #  ######    ##     #####  #    #  ######  #####
! #    #  #        #  #      #    #    #  #       #    #
! #    #  #####   #    #     #    ######  #####   #    #
! # ## #  #       ######     #    #    #  #       #####
! ##  ##  #       #    #     #    #    #  #       #   #
! #    #  ######  #    #     #    #    #  ######  #    #
!
! ---------------------------------------------------------------------
#ifdef SINGLE_POINT_MODEL
! This is the 0D version of the diurnal subroutine (formerly known as diurnalF)
! Copied on its entirety from Hewlley Imbuzeiro's 0D version. -fzm
subroutine diurnal (time, jday, test, dimforc)
      use inland_parameters, only: grav, rair, dtime, nband, pi
      use inland_comatm, only: snowa, raina, solad, solai, ta, psurf, qa, fira,&
                               rh, cloud, ua, coszen
      use inland_comforc, only: xta, xqa, xcld, xua, xprec, xlati, xsin, xlati,&
                                xlin
      use inland_comsatp

      implicit none

! Arguments
      integer jday, dimforc, time
      real*8 test

! determine the length of a precipitation event (between 4 and 24 hours),
! and time to start and end precipitation event. plen is in timesteps, while
! plens, startp, and endp are in seconds
 
! local variables
      integer ib,       & ! waveband number 1=visible, 2 near-IR
              it0

! FIXME: 'frac' variable name is already used in inland_comveg. This can turn 
!       out to be very misleading! -fzm
      real*8  rtime,    & ! time in hours
              orbit,    & ! earth's orbital angle (around the sun) in radians
              angle,    & ! solar hour angle in radians
              xdecl,    & ! solar declination angle
              sw,       & ! effective solar constant
              trans,    & ! solar transmission through the atmosphere
              fdiffuse, & ! fraction of indirect (diffuse) solar radiation
              frac,     & ! fraction of energy in each waveband
              dtair, rwork, pa, tdew, emb, ea, ec ! Hewlley did not say what
                                                  ! these are for. - fzm
#define WGEN_COMSAT
#include "inland_comsat.h"

! define working variables
      rwork = (grav / rair / 0.0065)

! * * *  ta, qa, rh, ua, cloud calculations * * *
      it0 = dimforc - test + 1
      ta(1) = xta(1,it0) + 273.15
      rh(1) = xqa(1,it0)
      cloud(1) = xcld(1,it0) / 100
      ua(1) = xua(1,it0)

! * * *  raina and snowa calculations * * *
! reset snow and rain to zero
      snowa(1) = 0.
      raina(1) = 0.
 
! Calcul of the actual raina, and snowa
      if (xta(1,it0) .gt. 2.5) then
         raina(1) = xprec(1,it0) / dtime
      else
         snowa(1) = xprec(1,it0) / dtime
      endif

! * * * calendar and orbital calculations * * *
! calculate time in hours
      rtime = (time - dtime / 2) / 3600.0

! calculate the earth's orbital angle (around the sun) in radians
      orbit = 2.0 * pi * float(jday) / 365.2425

! calculate the solar hour angle in radians
      angle  = 2.0 * pi * (rtime - 12.0) / 24.0
 
! calculate the current solar declination angle
! ref: global physical climatology, hartmann, appendix a
      xdecl = 0.006918 - 0.399912 * cos(orbit) + 0.070257 * sin(orbit) -  &
              0.006758 * cos(2.0 * orbit) + 0.000907 * sin(2.0 * orbit) - &
              0.002697 * cos(3.0 * orbit) + 0.001480 * sin(3.0 * orbit)
 
! calculate the effective solar constant, including effects of eccentricity
! ref: global physical climatology, hartmann, appendix a
      sw = 1370. * (1.000110 + 0.034221 * cos(orbit) + 0.001280 * sin(orbit) + &
           0.000719 * cos(2.0 * orbit) + 0.000077 * sin(2.0 * orbit)) 

! * * * solar calculations * * *
! calculate the cosine of the solar zenith angle
      coszen(1) = max(dble(0.0), (sin(xlati(1)) * sin(xdecl) + &
                  cos(xlati(1)) * cos(xdecl) * cos(angle)))

! calculate the solar transmission through the atmosphere
! using simple linear function of tranmission and cloud cover
!
! note that the 'cloud cover' data is typically obtained from
! sunshine hours -- not direct cloud observations
!
! where, cloud cover = 1 - sunshine fraction 
!
! different authors present different values for the slope and 
! intercept terms of this equation
!
! Friend, A: Parameterization of a global daily weather generator for
! terrestrial ecosystem and biogeochemical modelling, Ecological 
! Modelling
!
! Spitters et al., 1986: Separating the diffuse and direct component
! of global radiation and its implications for modeling canopy
! photosynthesis, Part I: Components of incoming radiation,
! Agricultural and Forest Meteorology, 38, 217-229.
!
! A. Friend       : trans = 0.251 + 0.509 * (1.0 - cloud(i))
! Spitters et al. : trans = 0.200 + 0.560 * (1.0 - cloud(i))
!
! we are using the values from A. Friend
      trans = 0.251 + 0.509 * (1.0 - cloud(1))

! calculate the fraction of indirect (diffuse) solar radiation
! based upon the cloud cover
!
! note that these relationships typically are measured for either
! monthly or daily timescales, and may not be exactly appropriate
! for hourly calculations -- however, in inland, cloud cover is fixed
! through the entire day so it may not make much difference
!
! method i --
!
! we use a simple empirical relationships from Nikolov and Zeller (1992)
!
! Nikolov, N. and K.F. Zeller, 1992:  A solar radiation algorithm for ecosystem
! dynamics models, Ecological Modelling, 61, 149-168.
      fdiffuse = 1.0045 + 0.0435 * trans - 3.5227 * trans**2 + 2.6313 * trans**3
 
      if (trans.gt.0.75) fdiffuse = 0.166
 
! method ii --
!
! another method was suggested by Spitters et al. (1986) based on
! long-term data from the Netherlands
!
! Spitters et al., 1986: Separating the diffuse and direct component
! of global radiation and its implications for modeling canopy
! photosynthesis, Part I: Components of incoming radiation,
! Agricultural and Forest Meteorology, 38, 217-229.
!       if ((trans.eq.0.00).and.(trans.lt.0.07)) then
!         fdiffuse = 1.0
!       else if ((trans.ge.0.07).and.(trans.lt.0.35)) then
!         fdiffuse = 1.0 - 2.3 * (trans - 0.07)**2
!       else if ((trans.ge.0.35).and.(trans.lt.0.75)) then
!         fdiffuse = 1.33 - 1.46 * trans
!       else
!         fdiffuse = 0.23
!       endif
!
! do for each waveband
      do 120 ib = 1, nband
 
! calculate the fraction in each waveband
!        frac = 0.46 + 0.08 * float(ib - 1)

! FIXME: there is absolutely no explanation on why are these changed from the
!       original version and hardcoded on the model.
         frac = 0.427 + 0.146 * float(ib - 1)

! calculate the direct and indirect solar radiation
!        solad(1,ib) = sw * coszen(1) * frac * trans * (1. - fdiffuse)  
!        solai(1,ib) = sw * coszen(1) * frac * trans * fdiffuse
 
         solad(1,ib) = xsin(1,it0) * frac * (1. - fdiffuse)
         solai(1,ib) = xsin(1,it0) * frac * fdiffuse
120   continue

! * * *  qa calculations * * *
      pa = 99208.7369088147 - (21.7473169490388 * xta(1,it0))

      qa(1) = rh(1)/100 * qsat(esat(ta(1)), pa)

! * * * ir flux calculations * * *
! clear-sky emissivity as a function of water vapor pressure
! and atmospheric temperature

! calculate the ir emissivity of the clear sky
! using equation from idso (1981) water resources res., 17, 295-304
      emb = 0.01 * (psurf(1) * qa(1) / (0.622 + qa(1)))
      ea  = 0.70 + 5.95e-5 * emb * exp (1500.0 / ta(1))

! assign the ir emissivity of clouds (assume they are ~black in the ir)
      ec = 0.950

! assign the temperature difference of emissions (air + cloud) from
! the surface air temperature
      if (coszen(1).le.0.) then
         dtair   = -8.0
      else
         dtair   = 15. * sin((rtime - 6) * (pi / 12))
      end if

! calcul of the dew point temperature in Celcius at 64 m
      tdew = (-237.3 * dlog(emb / 6.1078)) / (dlog(emb / 6.1078) - 17.269)

! total downward ir is equal to the sum of:
!
! (1) clear sky contribution to downward ir radiation flux
! (2) cloud contribution to downward ir radiation flux
!
!      fira(1) = (1. -  cloud(1)) * ea * stef * (ta(1) - dtair)**4 +
!     >                   cloud(1)  * ec * stef * (tdew + 273.16)**4 +
!     >                   12.
      fira(1) = xlin(1,it0)

#else /* SINGLE_POINT_MODEL */
! This is the ordinary 'diurnal' subroutine for the 2D weather generator.
      subroutine diurnal (time, jday, plens, startp, endp, starti, &
                          endi, seed,seed2,seed3,seed4)
! ---------------------------------------------------------------------
! common modules
      use inland_parameters
      use inland_comatm
      use inland_comwork
      use inland_comsatp
      use inland_subgrid
      use inland_comcrop
      use inland_comnitr
      use inland_combcs
      use inland_comveg
      use inland_control, only:imonth
 
      implicit none
 
! Arguments
      integer :: jday, seed, seed2, seed3(32), seed4, time ! current day
                
      real*8 :: plens,    & ! length of the precip event (s)
                startp,   & ! time to start precip event (s)
                endp,     & ! time to end precip event (s)
                starti,   &
                endi
 
! determine the length of a precipitation event (between 4 and 24 hours),
! and time to start and end precipitation event. plen is in timesteps, while
! plens, startp, and endp are in seconds
 
! local variables
      integer :: i, jj,   & ! loop indice
                 ib,      & ! waveband number 1= visible, 2= near-IR
                 ilpt       ! mlpt index

      real*8 :: gamma, qmin, qmax, qsa, emb, ea, ec, dtair, dtcloud, &
                rtime,      & ! time in hours
                orbit,      & ! earth's orbital angle (around the sun) in radians
                angle,      & ! solar hour angle in radians
                xdecl,      & ! solar declination angle
                sw,         & ! effective solar constant
                xlat,       & ! latitude in radians
                trans,      & ! solar transmission through the atmosphere
                fdiffuse,   & ! fraction of indirect (diffuse) solar radiation
                fracw         ! fraction of energy in each waveband
 
      integer :: checkP, niter, plen, plenmin, plenmax

! external random number generator -- must be real single precision!
      real ran2
 
#define WGEN_COMSAT
#include "inland_comsat.h"

! ---------------------------------------------------------------------- 
! * * * calendar and orbital calculations * * *
! ---------------------------------------------------------------------- 
 
! calculate time in hours
      rtime = time / 3600.0
 
! calculate the earth's orbital angle (around the sun) in radians
      orbit = 2.0 * pi * float(jday) / 365.2425

! calculate the solar hour angle in radians
      angle  = 2.0 * pi * (rtime - 12.0) / 24.0

! calculate the current solar declination angle
! ref: global physical climatology, hartmann, appendix a
      xdecl = 0.006918 - 0.399912 * cos(orbit) + 0.070257 * sin(orbit) -  &
              0.006758 * cos(2.0 * orbit) + 0.000907 * sin(2.0 * orbit) - &
              0.002697 * cos(3.0 * orbit) + 0.001480 * sin(3.0 * orbit)
 
! calculate the effective solar constant, including effects of eccentricity
! ref: global physical climatology, hartmann, appendix a
      sw = 1370. * (1.000110 + 0.034221 * cos(orbit) + 0.001280 * sin(orbit) + &
           0.000719 * cos(2.0 * orbit) + 0.000077 * sin(2.0 * orbit)) 
 
!9001  continue ! called by a 'goto'

! do for all gridcells
      do 100 i = 1, npoi
! ---------------------------------------------------------------------- 
! * * * solar calculations * * *
! ---------------------------------------------------------------------- 
! calculate the latitude in radians
          jj = latindex(i)
          xlat = latscale(jj) * pi / 180.0
 
! calculate the cosine of the solar zenith angle
          coszen(i) = max (dble(0.0), (sin(xlat) * sin(xdecl) + cos(xlat) * &
                      cos(xdecl) * cos(angle)))

! find daylength to be used in pheno subroutine
!
        daylength(i) = (180./pi)*((2.*60.)/15.)*(acos((coszen(i) &
                     - (sin(xlat)*sin(xdecl))) / (cos(xlat)*cos(xdecl))))

! calculate the solar transmission through the atmosphere
! using simple linear function of tranmission and cloud cover
! 
! note that the 'cloud cover' data is typically obtained from
! sunshine hours -- not direct cloud observations
! 
! where, cloud cover = 1 - sunshine fraction 
! 
! different authors present different values for the slope and 
! intercept terms of this equation
! 
! Friend, A: Parameterization of a global daily weather generator for
! terrestrial ecosystem and biogeochemical modelling, Ecological 
! Modelling
! 
! Spitters et al., 1986: Separating the diffuse and direct component
! of global radiation and its implications for modeling canopy
! photosynthesis, Part I: Components of incoming radiation,
! Agricultural and Forest Meteorology, 38, 217-229.
! 
! A. Friend       : trans = 0.251 + 0.509 * (1.0 - cloud(i))
! Spitters et al. : trans = 0.200 + 0.560 * (1.0 - cloud(i))
! 
! we are using the values from A. Friend
         trans = 0.251 + 0.509 * (1.0 - cloud(i)) 
 
! calculate the fraction of indirect (diffuse) solar radiation
! based upon the cloud cover
! 
! note that these relationships typically are measured for either
! monthly or daily timescales, and may not be exactly appropriate
! for hourly calculations -- however, in inland, cloud cover is fixed
! through the entire day so it may not make much difference
! 
! method i --
! 
! we use a simple empirical relationships from Nikolov and Zeller (1992)
! 
! Nikolov, N. and K.F. Zeller, 1992:  A solar radiation algorithm for ecosystem
! dynamics models, Ecological Modelling, 61, 149-168.
         fdiffuse = 1.0045 + 0.0435 * trans - 3.5227 * trans**2 + 2.6313 * &
                    trans**3
 
         if (trans.gt.0.75) fdiffuse = 0.166
 
! method ii --
! 
! another method was suggested by Spitters et al. (1986) based on
! long-term data from the Netherlands
! 
! Spitters et al., 1986: Separating the diffuse and direct component
! of global radiation and its implications for modeling canopy
! photosynthesis, Part I: Components of incoming radiation,
! Agricultural and Forest Meteorology, 38, 217-229.
! 
!       if ((trans.eq.0.00).and.(trans.lt.0.07)) then
!         fdiffuse = 1.0
!       else if ((trans.ge.0.07).and.(trans.lt.0.35)) then
!         fdiffuse = 1.0 - 2.3 * (trans - 0.07)**2
!       else if ((trans.ge.0.35).and.(trans.lt.0.75)) then
!         fdiffuse = 1.33 - 1.46 * trans
!       else
!         fdiffuse = 0.23
!       endif
! 
! do for each waveband
         do 120 ib = 1, nband
 
! calculate the fraction in each waveband
            fracw = 0.46 + 0.08 * float(ib - 1)
 
! calculate the direct and indirect solar radiation
            solad(i,ib) = sw * coszen(i) * fracw * trans * (1. - fdiffuse)  
            solai(i,ib) = sw * coszen(i) * fracw * trans * fdiffuse
120      continue
 
! ---------------------------------------------------------------------- 
! * * * temperature calculations * * *
! ---------------------------------------------------------------------- 
! 
! assign hourly temperatures using tmax and tmin 
! following Environmental Biophysics, by Campbell and Norman, p.23
! 
! this function fits a fourier series to the diurnal temperature cycle
! note that the maximum temperature occurs at 2:00 pm local solar time
! 
! note that the daily mean value of gamma is 0.44, 
! so td = 0.44 * tmax + 0.56 * tmin,  instead of
!    td = 0.50 * tmax + 0.50 * tmin
! 
         gamma = 0.44 - 0.46 * sin (pi / 12.0 * rtime + 0.9) + 0.11 * &
                 sin(2.0 * pi / 12.0 * rtime + 0.9)
         ta(i) = tmax(i) * gamma + tmin(i) * (1.0 - gamma)
 
! ---------------------------------------------------------------------- 
! * * * humidity calculations * * *
! ---------------------------------------------------------------------- 
! 
! adjust specific humidity against daily minimum temperatures
! 
! To do this, qa is written as an approximate sine function (same as ta)
! to preserve the daily mean specific humidity, while also preventing rh
! from exceeding 99% at night
! 
! Note that the daily mean RH is *not* preserved, and therefore the
! output RH will be slightly different from what was read in.
! 
! first adjust so that maximum RH cannot exceed 99% at night
         qmin = min (qd(i), 0.99 * qsat(esat(tmin(i)), psurf(i)))
         qmax = (qd(i) - 0.56 * qmin) / 0.44

! if needed, adjust again to 99% at other times of the day (in which
! case the daily mean *specific* humidity is also not preserved)
         qsa = 0.99 * qsat(esat(ta(i)), psurf(i))
 
! calculate the hourly specific humidity, using the above adjustments
         qa(i) = min (qsa, qmax * gamma + qmin * (1.0 - gamma))
 
! calculate the hourly relative humidity 
         rh(i) = 100.0 * qa(i) / qsat(esat(ta(i)), psurf(i))

! ---------------------------------------------------------------------- 
! * * * ir flux calculations * * *
! ---------------------------------------------------------------------- 
! 
! clear-sky emissivity as a function of water vapor pressure
! and atmospheric temperature
! 
! calculate the ir emissivity of the clear sky
! using equation from idso (1981) water resources res., 17, 295-304
         emb = 0.01 * (psurf(i) * qa(i) / (0.622 + qa(i)))
         ea  = 0.70 + 5.95e-5 * emb * exp (1500.0 / ta(i))
 
! assign the ir emissivity of clouds (assume they are ~black in the ir)
         ec = 0.950
 
! assign the temperature difference of emissions (air + cloud) from
! the surface air temperature
!
      if(isimagro .eq. 0) then
          dtair   = 0.0
          dtcloud = 0.0
      else
         dtair   = 2.0
         dtcloud = 2.0
      endif
!
! total downward ir is equal to the sum of:
! 
! (1) clear sky contribution to downward ir radiation flux
! (2) cloud contribution to downward ir radiation flux
         fira(i) = (1.-cloud(i)) * ea * stef * (ta(i)-dtair)**4 + cloud(i) * &
                   ec * stef * (ta(i) - dtcloud)**4

100   continue

! if mlpt > 1 do for all gridcells
! else
! do for non-tile gridcells because of multiple calls to ran2 will generate diff. results
! and then replicate ua, snowa, raina to subgrid tiles
! in any case, loop from 1 to npoi1, because when mlpt = 1, npoi1 = npoi
! TODO fix this when lbeg != 1 and lend != npoi1
      do 101 i = 1, npoi1
 
! ---------------------------------------------------------------------- 
! * * * wind speed calculations * * *
! ---------------------------------------------------------------------- 
! 
! following logic of the EPIC weather generator
! select random wind speed following this equation
         ua(i) = 1.13989 * ud(i) * (-log(ran2(seed,seed2,seed3,seed4)))**0.30 
 
! fix wind speeds to always be above 2.5 m/sec and below 10.0 m/sec
        if (isimagro .eq. 0) then
           ua(i) = max (dble(2.5), min (dble(10.0), ua(i)))
        else
           ua(i) = max (dble(0.2), min (dble(10.0), ua(i)))
        endif

9001  continue ! called by a 'goto'

! ---------------------------------------------------------------------- 
! * * * snow and rain calculations * * *
! ---------------------------------------------------------------------- 
 
! reset snow and rain to zero
         snowa(i) = 0.0
         raina(i) = 0.0
 
! determine the number of timesteps per day
         niter = int (86400.0 / dtime)
 
! change the rain length when the amount of rainfall/timestep is
! too high (at the first time step)
! FIXME ET - should this be done for every grid point? this doesn't make sense 

         if (time.lt.dtime) then
            plen = int( plens / dtime )
            plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
            plenmax = max (int (24.0 * 3600. / dtime), plenmin)
            checkP = 0

            do 85 while (((precip(i)/plen) .gt. 15).and.(plen.lt.plenmax))
               plen = plen + 1
               checkP = 1
85          continue

            if (checkP.eq.1) then
!              print *, 'WARNING: plen changed', i,
!     $             int(precip(i)), int(plens/dtime), plen
               plens = dtime * plen
               startp = dtime * min (niter-plen, int(ran2(seed,seed2,seed3,seed4)*(niter-plen+1)))
               endp = startp + plen *dtime
               goto 9001
            end if
         end if
 
! if precipitation event then calculate
         if (time.ge.startp .and. time.lt.endp) then  

! for rain / snow partitioning, make it all rain if 
! ta > 2.5 C and all snow if ta <= 2.5 C
! 
! reference:
! 
! Auer, A. H., 1974: The rain versus snow threshold temperatures,
! Weatherwise, 27, 67.
            if (ta(i)-273.15 .gt. 2.5) then
               raina(i) = precip(i) / plens
            else
               snowa(i) = precip(i) / plens
            endif
         endif
101   continue

! if using subgrid tiles, replicate ua, snowa, raina to subgrid tiles
if ( mlpt .gt. 1 ) then
      do i = 1, npoi1  
         do ilpt = 2,mlpt
            jj = subgrid_get_index(i,ilpt)
            if ( jj .ne. 0 ) then
               ua(jj) = ua(i)
               snowa(jj) = snowa(i)
               raina(jj) = raina(i)
            end if
         end do
      end do
endif

#endif /* SINGLE_POINT_MODEL */
      return
end subroutine diurnal
