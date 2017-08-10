#include "inland_config.h"
! ---------------------------------------------------------------------------------
subroutine vernalization(jday)
! ---------------------------------------------------------------------------------
!
! * * * only call subroutine for winter wheat * * *
!
! subroutine calculates vernalization and photoperiod effects on gdd accumulation
! in winter wheat varieties. Thermal time accumulation is reduced in 1st period until
! plant is fully vernalized. During this time of emergence to spikelet formation,  
! photoperiod can also have a drastic effect on plant development.
!
!
!
! common blocks
!
      use inland_parameters
      use inland_control,only: iyear,imonth,iday
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
      use inland_comsno
!
! local variables
!
      real*8 p1d, p1v,    &
             tcrown,      &
             vd,vd1,vd2,  &
             tkil,tbase,  &
             hti 
!
      integer  i,       &
               jday
 
! photoperiod factor calculation      
! genetic constant - can be modified
!
      p1d = 0.004  ! average for genotypes from Ritchey, 1991.
!                    ! Modeling plant and soil systems: Wheat phasic development
      p1v = 0.003  ! average for genotypes from Ritchey, 1991.  
!
      do 100 i = lbeg, lend
!
! only go through routine if winter wheat has been planted, is living,
! and the factor is not already 1.0
! 


         if (croplive(i,15) .eq. 1.0 .and. vf(i) .ne. 1.0) then
!
! from CERES-wheat J.T. Ritchey
!
!          jj   = latindex(i)
!          xlat = latscale(jj) 
!          s1 = sin(xlat * 0.01745)
!          c1 = cos(xlat * 0.01745)
!          dec  = 0.4093 * sin(0.0172 * (jday - 82.2)) 
!          dlv  = ((-s1 * sin(dec) - 0.1047)/(c1 * cos(dec)))

!          if (dlv .lt. -0.87) dlv = -0.87
!             hrlt = 7.639 * acos(dlv)  
!             df   = 1 - p1d * (20.0 - hrlt)**2
!
! daylength calculation in IBIS is in minutes - convert to hours
!
            df(i) = 1 - p1d * (20.0 - daylength(i)/60.0)**2

!
! for all equations - temperatures must be in degrees (C)
! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
! snow depth in centimeters
!
            if (td(i) .lt. 273.16) then
               tcrown = 2.0 + (td(i) - 273.16) * (0.4 + 0.0018 * &
                       (min(adsnod(i)*100., 15.0) - 15.0)**2)
            else
               tcrown = td(i) - 273.16 
            endif
! 
! vernalization factor calculation
! if vf(i) = 1.  then plant is fully vernalized - and thermal time
! accumulation in phase 1 will be unaffected 
! cumulative vd will need to be reset to 0 at planting of crop
!
            if (vf(i) .ne. 1.0) vd = 0.
            if (vf(i) .ne. 1.0 .and. tmax(i) .gt. 273.16) then
               if (tmin(i) .le. 288.16) then
                  vd1 = 1.4 - 0.0778 * (tcrown)
                  vd2 = 0.5 + 13.44/((tmax(i) - tmin(i) + 3.0)**2) * tcrown
                  vd = max(0.0, min(1., vd1, vd2))
                  cumvd(i) = cumvd(i) + vd
               endif
!
               if (cumvd(i) .lt. 10 .and. tmax(i) .gt. 303.16) then
                  cumvd(i) = cumvd(i) - 0.5 * (tmax(i) - 303.16)   
               endif 
               cumvd(i) = max(0.0, cumvd(i))    ! must be greater than 0 
               vf(i) = 1.0 - p1v * (50.0 - cumvd(i))
               vf(i) = max(0.0, min(vf(i), 1.0))       ! must be between 0 - 1
            endif
!             
! calculate cold hardening of plant
! determines for winter wheat varieties whether the plant has completed
! a period of cold hardening to protect it from freezing temperatures.  If  it has
! not, then exposure could result in death or killing of plants.
!
! there are two distinct phases of hardening 
!
         tbase = 0.0
         if (tmin(i) .le. 270.16 .or. hdidx(i) .ne. 0) then
            hti =  1.0
            if (hdidx(i) .ge. hti) then   ! done with phase 1
               hdidx(i) = hdidx(i) + 0.083
               if (hdidx(i) .gt. hti * 2.0) hdidx(i) = hti * 2.0  
            endif
            if (tmax(i) .ge. tbase + 283.16) then
               hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16)
               if (hdidx(i) .gt. hti) hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16) 
                  hdidx(i) = max(0.0, hdidx(i))
               endif
            elseif (tcrown .ge. tbase - 1.0) then
               if (tcrown .le. tbase + 8.0) then
                  hdidx(i) = hdidx(i) + 0.1 - ((tcrown - tbase + 3.5)**2 / 506.0)
                  if (hdidx(i) .ge. hti .and. tcrown .le. tbase + 0.0) then
                     hdidx(i) = hdidx(i) + 0.083
                     if (hdidx(i) .gt. hti * 2.0) hdidx(i) = hti * 2.0
                  endif
               endif 
               if (tmax(i) .ge. tbase + 283.16) then
                  hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16)
                  if (hdidx(i) .gt. hti) hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16) 
                     hdidx(i) = max(0.0, hdidx(i))
               endif
            endif
!      
! calculate what the wheat killing temperature 
! there is a linear inverse relationship between 
! hardening of the plant and the killing temperature or
! threshold that the plant can withstand 
! when plant is fully-hardened (hdidx = 2), the killing threshold is -18 C
!
! will have to develop some type of relationship that reduces LAI and
! biomass pools in response to cold damaged crop 
!
            if (tmin(i) .le. 267.16) then
               tkil = (tbase - 6.0) - 6.0 * hdidx(i)         
               if (tkil .ge. tcrown) then
                  if ((0.95 - 0.02 * (tcrown - tkil)**2) .ge. 0.02) then
                     write (*,*)  'crop damaged by cold temperatures'
                  else
                     croplive(i,15) = 0.0   ! kill winter wheat
                     write (*,*)  '95% of crop killed by cold temperatures'
                  endif
               endif
            endif
            else if (croplive(i,13) .eq. 1 .or. croplive(i,14) .eq. 1) then
               vf(i) = 1.0
         endif  ! croplive, vf check
!
100   continue
!
      return
end subroutine vernalization

