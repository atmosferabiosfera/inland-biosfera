#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine sumyear(mcsec, loopi, kpti, kptj)
! ---------------------------------------------------------------------
      use inland_com1d
      use inland_parameters
      use inland_control, only: iday, imonth, isimveg, isimfire
      use inland_comatm
      use inland_comhyd
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comfire
      use inland_comnitr
      use inland_comcrop, only:isimagro

      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      integer loopi    ! index of little vector in big vector
      integer kpti     ! index of 1st point of little vector in big lpt vector
      integer kptj     ! index of last point of little vector

      integer mcsec    ! current second in day (passed in)

! local variables
      real*8 rwork, rwork2, rwork3, rwork4, & !
             rdepth,   &! 1/total soil depth over 4 1st layers
             solartot, &! total incoming radiation (direct + diffuse, visible + nearIR)
             reflectot,&! total incoming radiation (direct + diffuse, visible + nearIR)
             soiltemp, &! average soil temp for 4 1st layers
             soilmois, &! average soil moisture for 4 1st layers 
             soilice,  &! average soil ice for 4 1st layers 
             vwc,      &! total liquid + ice content of 4 1st layers
             awc,      &! total available water (+ ice) content of 4 1st layer
             water,    &! fire factor: total water content of 1st layer (liquid+ice)
             waterfrac,&! fire factor: available water content of 1st layer
             fueldry,  &! fire factor
             allroots   ! annual average root biomass

      integer i, k, np, nytimes

      logical lastts ! is this the last timestep of a given month?

! ---------------------------------------------------------------------
! * * * update counters and working variables * * *
! ---------------------------------------------------------------------
!
! reset sumyear if the first timestep of the year
      if ((mcsec.eq.0).and.(iday.eq.1).and.(imonth.eq.1)) then
         nytime(loopi) = 0

         ! TODO reset other vars as well
         ayburnfrac(:) = 0
         amburnpft(:,:) = 0
      end if

      lastts = .false.
      if ( ( mcsec .eq. (86400 - dtime) ) .and. (iday .eq. ndaypm(imonth) ) ) lastts = .true.
      
! accumulate yearly output
      nytime(loopi) = nytime(loopi) + 1

! working variables
      nytimes = nytime(loopi)

! rwork4 is for nitrogen mineralization conversion
      rwork  = 1. / float(nytimes)
      rwork2 = float(ndaypy) * 86400.
      rwork3 = float(ndaypy) * 86400. * 12.e-3
      rwork4 = float(ndaypy) * 86400. * 14.e-3
      rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))

! ---------------------------------------------------------------------
! begin global grid
! ---------------------------------------------------------------------
      do 100 i = kpti, kptj

! ---------------------------------------------------------------------
! * * * annual energy budget terms * * *
! ---------------------------------------------------------------------
         solartot = solad(i,1) + solad(i,2) + solai(i,1) + solai(i,2)
         reflectot = asurd(i,1) * solad(i,1) + asurd(i,2) * solad(i,2) + &
                     asuri(i,1) * solai(i,1) + asuri(i,2) * solai(i,2)
         aysolar(i) = ((nytimes-1) * aysolar(i) + solartot) * rwork
         ayreflect(i) = ((nytimes-1) * ayreflect(i) + reflectot) * rwork

! normalize the albedo calcuation
!        ayalbedo(i) = ayalbedo(i) / (aysolar(i) + epsilon)
         ayirup(i) = ((nytimes-1) * ayirup(i) + firb(i)) * rwork
         ayirdown(i) = ((nytimes-1) * ayirdown(i) + fira(i)) * rwork
         aysens(i) = ((nytimes-1) * aysens(i) - fsena(i)) * rwork
         aylatent(i) = ((nytimes-1) * aylatent(i) - fvapa(i) * hvap) * rwork

! ---------------------------------------------------------------------
! * * * annual water budget terms * * *
! ---------------------------------------------------------------------
        ayprcp(i)=((nytimes-1)*ayprcp(i) + (raina(i) + snowa(i))*rwork2)*rwork
        ayaet(i) = ((nytimes-1) * ayaet(i) - fvapa(i) * rwork2) * rwork
        aytrans(i) = ((nytimes-1) * aytrans(i) + gtrans(i) * rwork2) * rwork
        aysrunoff(i)=((nytimes-1) * aysrunoff(i) + grunof(i) * rwork2) * rwork
        aydrainage(i)=((nytimes-1)*aydrainage(i) + gdrain(i)*rwork2)*rwork
        if(isimagro .eq. 0) then
           aytrunoff(i) = aysrunoff(i) + aydrainage(i) 
        else
           aytrunoff(i)  = ((nytimes-1) * aytrunoff(i) + (grunof(i) + gdrain(i)) * rwork2) * rwork
        endif
! estimate the change in total surface water content (intercepted
! water and snow, soil water and ice, snow).
         wtotp(i) = wtot(i)
         wtot(i) = (wliqu(i)+wsnou(i)) * fu(i) * 2.0 * lai(i,2) + &
                   (wliqs(i)+wsnos(i)) * fu(i) * 2.0 * sai(i,2) + &
                   (wliql(i)+wsnol(i)) * fl(i) * 2.0 * &
                   (lai(i,1) + sai(i,1)) * (1. - fi(i))
         wtot(i) = wtot(i) + wpud(i) + wipud(i)
         do k = 1, nsoilay
            wtot(i) = wtot(i) + &
            poros(i,k)*wsoi(i,k) * (1.-wisoi(i,k)) * hsoi(k) * rhow + &
            poros(i,k)*wisoi(i,k) * hsoi(k)*rhow
         end do
         do k = 1, nsnolay
            wtot(i) = wtot(i) + fi(i) * rhos * hsno(i,k)
         end do
         dwtot(i) = ((nytimes-1) * dwtot(i) + wtot(i) - wtotp(i)) * rwork

! ---------------------------------------------------------------------
! * * * fire parameters * * *
! ---------------------------------------------------------------------
if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
         aypbio(i) = ((nytimes-1) * aypbio(i) + pbio(i)) * rwork
         aypmoi(i) = ((nytimes-1) * aypmoi(i) + pmoi(i)) * rwork
         aypign(i) = ((nytimes-1) * aypign(i) + Pign(i)) * rwork
         aypfire(i) = ((nytimes-1) * aypfire(i) + Pfire(i)) * rwork
         aysrate(i) = ((nytimes-1) * aysrate(i) + srate(i)) * rwork
         ayabday(i) = ((nytimes-1) * ayabday(i) + abday(i)) * rwork
         pfireyr(i) = aypfire(i)

         ! only update burnfrac at last timestep of the month
         !ayburnarea(i) = ((nytimes-1) * ayburnarea(i) + burnarea(i)) * rwork
         !ayburnpft(i,np) = ((nytimes-1) * ayburnpft(i,np) + burnpft(i,np)) * rwork
         if ( lastts ) then
            ayburnfrac(i) = ayburnfrac(i) + amburnfrac(i)
            do np = 1, npft
               ayburnpft(i,np) = ayburnpft(i,np) + amburnpft(i,np)
            end do
         end if

end if
! ---------------------------------------------------------------------
! * * * annual soil parameters * * *
! ---------------------------------------------------------------------
         soiltemp = 0.0
         soilmois = 0.0
         soilice  = 0.0
         vwc = 0.0
         awc = 0.0

! averages for first 4 layers of soil
         do 110 k = 1, 4
            soiltemp =  soiltemp + tsoi(i,k)  * hsoi(k)
            soilmois =  soilmois + wsoi(i,k)  * hsoi(k)
            soilice  =  soilice  + wisoi(i,k) * hsoi(k)
            vwc = vwc + (wisoi(i,k) + (1. - wisoi(i,k)) * wsoi(i,k)) * &
                  hsoi(k) * poros(i,k)
            awc = awc + max (dble(0.0), (wisoi(i,k) + &
                  (1. - wisoi(i,k)) * wsoi(i,k)) - swilt(i,k)) * &
                  hsoi(k) * poros(i,k) * 100.0
110      continue

! average soil and air temperatures
         soiltemp = soiltemp * rdepth - 273.16
         soilmois = soilmois * rdepth
         soilice  = soilice  * rdepth
         vwc = vwc * rdepth
         awc = awc * rdepth

! annual average soil moisture and soil ice
         aywsoi(i) = ((nytimes-1) * aywsoi(i) + soilmois) * rwork
         aywisoi(i) = ((nytimes-1) * aywisoi(i) + soilice) * rwork
         aytsoi(i) = ((nytimes-1) * aytsoi(i) + soiltemp) * rwork
         ayvwc(i) = ((nytimes-1) * ayvwc(i) + vwc) * rwork
         ayawc(i) = ((nytimes-1) * ayawc(i) + awc) * rwork

! soil moisture stress
       if(isimagro .eq. 0)then
!       aystresstu(i) = rwork * ((nytimes-1) * aystresstu(i) +
!    >                  stresstu(i))
!       aystresstl(i) = rwork * ((nytimes-1) * aystresstl(i) +
!    >                  stresstl(i))
       else
            aystresstu(i) = rwork * ((nytimes-1) * aystresstu(i) + &
                            stresstu(i))
            aystresstl(i) = rwork * ((nytimes-1) * aystresstl(i) + &
                            stresstl(i))
       endif

! ---------------------------------------------------------------------
! * * * determine annual gpp * * *
! ---------------------------------------------------------------------
!
! gross primary production of each plant type
         do np = 1, npft
            aygpp(i,np) = ((nytimes-1) * aygpp(i,np) + &
                          tgpp(i,np)  * rwork3) * rwork
         end do

! gross primary production of the entire gridcell
       if(isimagro .gt. 0) then
         aygpptot(i) = aygpp(i,1) + aygpp(i,2) + aygpp(i,3)    + &
                       aygpp(i,4) + aygpp(i,5) + aygpp(i,6)    + &
                       aygpp(i,7) + aygpp(i,8) + aygpp(i,9)    + &
                       aygpp(i,10) + aygpp(i,11) + aygpp(i,12) + &
                       aygpp(i,13)+ aygpp(i,14) + aygpp(i,15)  + &
                       aygpp(i,16)
       else
         aygpptot(i) = aygpp(i,1) + aygpp(i,2) + aygpp(i,3) + &
                       aygpp(i,4) + aygpp(i,5) + aygpp(i,6) + &
                       aygpp(i,7) + aygpp(i,8) + aygpp(i,9) + &
                       aygpp(i,10) + aygpp(i,11) + aygpp(i,12)

       endif
! ---------------------------------------------------------------------
! * * * determine annual npp * * *
! ---------------------------------------------------------------------
!
! net primary production of each plant type
         do np = 1, npft
            aynpp(i,np) = ((nytimes-1) * aynpp(i,np) + tnpp(i,np) * rwork3) * &
                          rwork
         end do

! net primary production of the entire gridcell
       if(isimagro .gt. 0) then
         aynpptot(i) = aynpp(i,1) + aynpp(i,2) + aynpp(i,3)    + &
                       aynpp(i,4) + aynpp(i,5) + aynpp(i,6)    + &
                       aynpp(i,7) + aynpp(i,8) + aynpp(i,9)    + &
                       aynpp(i,10) + aynpp(i,11) + aynpp(i,12) + &
                       aynpp(i,13) + aynpp(i,14) + aynpp(i,15) + &
                       aynpp(i,16)
       else
         aynpptot(i) = aynpp(i,1) + aynpp(i,2) + aynpp(i,3) + &
                       aynpp(i,4) + aynpp(i,5) + aynpp(i,6) + &
                       aynpp(i,7) + aynpp(i,8) + aynpp(i,9) + &
                       aynpp(i,10) + aynpp(i,11) + aynpp(i,12)
       endif

! ---------------------------------------------------------------------
! * * * annual carbon budget terms * * *
! ---------------------------------------------------------------------
!
! fire factor used in vegetation dynamics calculations
         water     = wisoi(i,1) + (1. - wisoi(i,1)) * wsoi(i,1)
         waterfrac = (water - swilt(i,1)) / (1. - swilt(i,1))
         fueldry = max (dble(0.0), min (dble(1.0), -2.0 * (waterfrac - 0.5)))
         firefac(i) = ((nytimes-1) * firefac(i) + fueldry) * rwork

! increment annual total co2 respiration from microbes
! tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
         ayco2mic(i) = ((nytimes-1) * ayco2mic(i) + tco2mic(i) * rwork3) * rwork

! increment annual total co2 respiration from roots
         ayco2root(i)=((nytimes-1)*ayco2root(i) + tco2root(i)*rwork3)*rwork

! calculate annual total co2 respiration from soil
         ayco2soi(i) = ayco2root(i) + ayco2mic(i)

!  annual net ecosystem co2 flux - npp total minus microbial respiration
!  the npp total includes losses from root respiration
         ayneetot(i) = aynpptot(i) - ayco2mic(i)


! annual average root biomass
       if(isimagro .gt. 0) then
         allroots = cbior(i,1) + cbior(i,2)  + cbior(i,3)  + cbior(i,4)  + &
                    cbior(i,5) + cbior(i,6)  + cbior(i,7)  + cbior(i,8)  + &
                    cbior(i,9) + cbior(i,10) + cbior(i,11) + cbior(i,12) + &
                    cbior(i,13) + cbior(i,14) + cbior(i,15) + cbior(i,16)
       else
         allroots = cbior(i,1) + cbior(i,2)  + cbior(i,3)  + cbior(i,4) + &
                    cbior(i,5) + cbior(i,6)  + cbior(i,7)  + cbior(i,8) + &
                    cbior(i,9) + cbior(i,10) + cbior(i,11) + cbior(i,12)
       endif
         ayrootbio(i) = ((nytimes-1) * ayrootbio(i) + allroots) * rwork

! ---------------------------------------------------------------------
! * * * annual biogeochemistry terms * * *
! ---------------------------------------------------------------------
!
! increment annual total of net nitrogen mineralization
! value for tnmin is calculated in biogeochem.f
         aynmintot(i) = ((nytimes-1) * aynmintot(i) + tnmin(i) * rwork4) * rwork

     if(isimagro .gt. 0) then
!
! mineralization (gross)
!
        aymintot(i)  = ((nytimes-1) * aymintot(i) + &
                     totmin(i) * rwork4) * rwork
!
! immobilization (gross)
!
        ayimmtot(i)  = ((nytimes-1) * ayimmtot(i) + &
                     totimm(i) * rwork4) * rwork 
!
! other mineralization/immobilization
! from non-microbial transformations
!
        aynreltot(i) = ((nytimes-1) * aynreltot(i) + &
                     totnrel(i) * rwork4) * rwork 
!
     endif ! check for crop existence
! other biogeochemistry variables
         ayalit(i)  = ((nytimes-1) * ayalit(i)  + totalit(i))  * rwork
         ayblit(i)  = ((nytimes-1) * ayblit(i)  + totrlit(i))  * rwork
         aycsoi(i)  = ((nytimes-1) * aycsoi(i)  + totcsoi(i))  * rwork
         aycmic(i)  = ((nytimes-1) * aycmic(i)  + totcmic(i))  * rwork
         ayanlit(i) = ((nytimes-1) * ayanlit(i) + totanlit(i)) * rwork
         aybnlit(i) = ((nytimes-1) * aybnlit(i) + totrnlit(i)) * rwork
         aynsoi(i)  = ((nytimes-1) * aynsoi(i)  + totnsoi(i))  * rwork
100   continue

      return
end subroutine sumyear
