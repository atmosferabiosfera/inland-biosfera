#include "inland_config.h"
! ---------------------------------------------------------------------
      subroutine diurnalmet (timed, jday, plens, startp, endp, &
                 irrigate, ilens, starti, endi)
! ---------------------------------------------------------------------
!
! common blocks
!
!
      use inland_comveg
      use inland_comhour
      use inland_comcrop
      use inland_comnitr
      use inland_parameters
      use inland_comatm
      use inland_comwork
      use inland_comsatp
 
      implicit none
 
! Arguments
!
      integer jday,    & ! current day
              irrigate, &
              timed
!
      real*8 plens,      & ! length of the precip event (s)
           startp,     & ! time to start precip event (s)
           endp,       & ! time to end precip event (s)
           ilens,      &
           starti,     &
           endi,       &
           truecloud 
!
      integer i,        & ! loop indice
              jj,       &
              ib          ! waveband number 1= visible, 2= near-IR
!
      real rtime,        & ! time in hours
           orbit,        & ! earth's orbital angle (around the sun) in radians
           angle,        & ! solar hour angle in radians
           xdecl,        & ! solar declination angle
           sw,           & ! effective solar constant
           xlat,         & ! latitude in radians
           trans,        & ! solar transmission through the atmosphere
           fdiffuse,     & ! fraction of indirect (diffuse) solar radiation
           fracw,        & ! fraction of energy in each waveband
           wfrac,        &
           gamma,        & !
           qmin,         & !
           qmax,         &
           qsa,          &
           ran2,         &
           emb,          &
           ea,           &
           ec,           &
           dtair,        &
           dtcloud
!
      integer checkP,    &
              niter,     &
              plen,      &
              plenmin,   &
              plenmax
!
#define WGEN_COMSAT
#include "inland_comsat.h"
!
!
! ---------------------------------------------------------------------- 
! * * * calendar and orbital calculations * * *
! ---------------------------------------------------------------------- 
!
! calculate time in hours
!

      rtime = timed / 3600.0
!
! calculate the earth's orbital angle (around the sun) in radians
!
      orbit = 2.0 * pi * float(jday) / 365.2425
!
! calculate the solar hour angle in radians
!
      angle  = 2.0 * pi * (rtime - 12.0) / 24.0
!
! calculate the current solar declination angle
! ref: global physical climatology, hartmann, appendix a
!
      xdecl =  0.006918                    -    &
               0.399912 * cos(orbit)       +    &
               0.070257 * sin(orbit)       -    &
               0.006758 * cos(2.0 * orbit) +    &
               0.000907 * sin(2.0 * orbit) -    &
               0.002697 * cos(3.0 * orbit) +    &
               0.001480 * sin(3.0 * orbit)
!
! calculate the effective solar constant, including effects of eccentricity
! ref: global physical climatology, hartmann, appendix a
!
      sw = 1370. * (1.000110 +                      &
                    0.034221 * cos(orbit)        +  &
                    0.001280 * sin(orbit)        +  &
                    0.000719 * cos(2.0 * orbit)  +  &
                    0.000077 * sin(2.0 * orbit)) 
!
! 9001 continue
!
! do for all gridcells
!

      do 100 i = lbeg, lend

!
! ---------------------------------------------------------------------- 
! * * * solar calculations * * *
! ---------------------------------------------------------------------- 
!
! calculate the latitude in radians
!
        jj = latindex(i)
!
        xlat = latscale(jj) * pi / 180.0
!
! calculate the cosine of the solar zenith angle
!
        coszen(i) = max (0.0, (sin(xlat) * sin(xdecl) +  &
                    cos(xlat) * cos(xdecl) * cos(angle)))
!
! find daylength to be used in pheno subroutine
!
        daylength(i) = (180./pi)*((2.*60.)/15.)*(acos((coszen(i) &
           - (sin(xlat)*sin(xdecl))) / (cos(xlat)*cos(xdecl))))
!
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
!
!cjk     trans = 0.251 + 0.509 * (1.0 - cloud(i)) 
!        trans = cloud(i) / sw     ! cloud(i) is surface insolation
        trans = cloud(i) / (sw*coszen(i)) 
	trans = max(0.,min(1.,trans))

!      if(i.eq.1)print*,jday,time/3600.,cloud(i),sw*coszen(i),trans
!
!cjk
!        if (jday .eq. 1) open (26, file='met.trans',status='unknown')
!        write(26,*) jday, istep, trans
!
! calculate the fraction of indirect (diffuse) solar radiation
! based upon the cloud cover
!
! note that these relationships typically are measured for either
! monthly or daily timescales, and may not be exactly appropriate
! for hourly calculations -- however, in ibis, cloud cover is fixed
! through the entire day so it may not make much difference
!
! method i --
!
! we use a simple empirical relationships from Nikolov and Zeller (1992)
!
! Nikolov, N. and K.F. Zeller, 1992:  A solar radiation algorithm for ecosystem
! dynamics models, Ecological Modelling, 61, 149-168.
!
        fdiffuse = 1.0045 + 0.0435 * trans &
                 - 3.5227 * trans**2       &
                 + 2.6313 * trans**3

!
        if (trans.gt.0.75) fdiffuse = 0.166
!
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
!
        do 120 ib = 1, nband
!
! calculate the fraction in each waveband
!
          wfrac = 0.46 + 0.08 * float(ib - 1)  !visible 0.46 and NIR 0.54
!
! calculate the direct and indirect solar radiation
!
! cjk       solad(i,ib) = sw * coszen(i) * wfrac * trans *
! cjk  >                  (1. - fdiffuse)  
! cjk
          solad(i,ib) = wfrac * cloud(i) * (1. - fdiffuse)
!
          solai(i,ib) = wfrac * cloud(i) * fdiffuse

!
  120   continue
!
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
        gamma = 0.44 - 0.46 * sin ( pi / 12.0 * rtime + 0.9) &
              + 0.11 * sin (2.0 * pi / 12.0 * rtime + 0.9)


!      if(time.eq.13*dtime.and.i.eq.1)print*,jday,time,ta(i),ta(i)
!
! cjk    ta(i) = tmax(i) * gamma + tmin(i) * (1.0 - gamma)
!
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
!
! cjk     qmin = min (qd(i), 0.99 * qsat(esat(tmin(i)), psurf(i)))
! cjk     qmax = (qd(i) - 0.56 * qmin) / 0.44
!
! if needed, adjust again to 99% at other times of the day (in which
! case the daily mean *specific* humidity is also not preserved)
!
        qsa  = 0.99 * qsat (esat (ta(i)), psurf(i))
!
! calculate the hourly specific humidity, using the above adjustments
!
! cjk     qa(i) = min (qsa, qmax * gamma + qmin * (1.0 - gamma))
        qa(i) = min (qsa, qd(i))
!
! calculate the hourly relative humidity 
!
        rh(i) = 100.0 * qa(i) / qsat (esat (ta(i)), psurf(i))

!
! ---------------------------------------------------------------------- 
! * * * wind speed calculations * * *
! ---------------------------------------------------------------------- 
!
! following logic of the EPIC weather generator
! select random wind speed following this equation
!
! cjk     ua(i) = 1.13989 * ud(i) * (-log(ran2(seed)))**0.30 
!
! fix wind speeds to always be above 2.5 m/sec and below 10.0 m/sec
!
! ---------------------------------------------------------------------- 
! * * * ir flux calculations * * *
! ---------------------------------------------------------------------- 
!
! clear-sky emissivity as a function of water vapor pressure
! and atmospheric temperature
!
! calculate the ir emissivity of the clear sky
! using equation from idso (1981) water resources res., 17, 295-304
!
        emb = 0.01 * (psurf(i) * qa(i) / (0.622 + qa(i)))
        ea  = 0.70 + 5.95e-5 * emb * exp (1500.0 / ta(i))

!      if(time.eq.13*dtime.and.i.eq.1)print*,jday,time,psurf(i),'  1'
!
! assign the ir emissivity of clouds (assume they are ~black in the ir)
!
        ec = 0.950
!
! assign the temperature difference of emissions (air + cloud) from
! the surface air temperature
!
        dtair   = 2.0
        dtcloud = 2.0
!
! total downward ir is equal to the sum of:
!
! (1) clear sky contribution to downward ir radiation flux
! (2) cloud contribution to downward ir radiation flux
!
! cjk 
!
        truecloud = 1. - ((trans - 0.251) / 0.509) 
        fira(i) = (1. -  truecloud) * ea * stef * (ta(i) - dtair  )**4 + &
                truecloud  * ec * stef * (ta(i) - dtcloud)**4

!	if(time.eq.13*dtime.and.i.eq.1)
!     >	print*,jday,time,fira(i)
!
!      fira(i) = (1. -  cloud(i)) * ea * stef * (ta(i) - dtair  )**4 +
!   >                   cloud(i)  * ec * stef * (ta(i) - dtcloud)**4
!
! ---------------------------------------------------------------------- 
! * * * snow and rain calculations * * *
! ---------------------------------------------------------------------- 
!
! reset snow and rain to zero
!
        snowa(i) = 0.0
        raina(i) = 0.0
!
! determine the number of timesteps per day
!
        niter = int (86400.0 / dtime)
!
! change the rain length when the amount of rainfall/timestep is
! too high (at the first time step)
!
!        if (time.lt.dtime) then
!
!           plen = plens / dtime
!           plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
!           plenmax = max (int (24.0 * 3600. / dtime), plenmin)
!           checkP = 0
!
!           do  while (((precip(i)/plen) .gt. 15).and.(plen.lt.plenmax))
!              plen = plen + 1
!              checkP = 1
!           end do
!
!           if (checkP.eq.1) then
!
!              print *, 'WARNING: plen changed', i,
!     $             int(precip(i)), int(plens/dtime), plen
!              plens = dtime * plen
!              startp = dtime * min (niter-plen,
!     >             int(ran2(seed)*(niter-plen+1)))
!              endp = startp + plen *dtime
!              goto 9001
!           end if
!
!        end if
!
! if precipitation event then calculate
!
! cjk        if (time.ge.startp .and. time.lt.endp) then  
!
! for rain / snow partitioning, make it all rain if 
! ta > 2.5 C and all snow if ta <= 2.5 C
!
! reference:
!
! Auer, A. H., 1974: The rain versus snow threshold temperatures,
! Weatherwise, 27, 67.
!
!

          if (ta(i)-273.15 .gt. 2.5) then
            raina(i) = precip(i) / plens
          else
            snowa(i) = precip(i) / plens
          endif


!      if(i.eq.1)print*,jday,time/3600.,precip(i),raina(i),plens
!
! cjk        endif
!
! ---------------------------------------------------------------------- 
! * * * irrigation calculations * * *
! ---------------------------------------------------------------------- 
!
! reset rate of irrigation application per timestep 
!
        xirriga(i) = 0.0
!
! if precipitation event - then no irrigation that day 
!
!
        if (timed.ge.starti .and. timed.lt.endi &
           .and. irrigate .eq. 1              &
           .and. precip(i) .eq. 0.00) then  
!
          xirriga(i) = xirrig(i) / ilens
!
! update annual total - totirrig
! rate of irrigation multiplied by length of timestep (mm/s * s) = total applied
! for this timestep 
!
          totirrig(i) = totirrig(i) + (xirriga(i) * dtime)
!
        endif
!
  100 continue
!
   return
end subroutine diurnalmet
