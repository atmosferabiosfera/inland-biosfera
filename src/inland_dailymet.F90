#include "inland_config.h"
! -----------------------------------------------------------------------------------
      subroutine dailymet (imonth, jday)
! -----------------------------------------------------------------------------------
!
! overview
!
! this routine generates daily weather conditions from hourly obsarvations -
!
!
! specifically, this routine needs daily values of
!
!  - daily total precipitation
!  - daily maximum temperature
!  - daily minimum temperature
!  - daily average cloud cover
!  - daily average relative humidity
!  - daily average wind speed
!
!
! common blocks
!
!
      use inland_parameters
!      use inland_control
      use inland_comveg
      use inland_combcs
      use inland_comatm
      use inland_comsoi
      use inland_comcrop
      use inland_comsatp
      use inland_comsum

      implicit none
!
! Arguments
!
      integer jday,    &     ! 1 if reading in daily weather data
              iday,    &
              imonth,  &
              iyear
!
! local variables
!
      integer it1w,    &   ! indice of previous month (interpolation)
              it2w,    &   ! indice of following month (interpolation)
              i,       &
              j,       &
              k            ! loop indice
!
      real rwork,      &   !
           omcloud,    &   ! cloud cover
           omqd,       &   ! humidity
           omtmax,     &   ! maximum temperature
           ran2,       &   ! function random number generator
           dt,         &   ! used for interpolation
           pwet,       &   ! monthly-average probability of rainy day
           pwd,        &   ! probability of a wet day after a dry day
           pww,        &   ! probability of a wet day after a wet day
           rndnum,     &   ! random number to decide if wet or dry day
           rainpwd,    &   ! average rainfall per wet day
           alpha,      &   ! parameter for gamma function
           beta,       &   ! parameter for gamma function
           aa,         &
           ab,         &
           tr1,        &
           tr2,        &
           rn1,        &   ! random number
           rn2,        &   ! random number
           rn3,        &   ! random number
           rn,         &   ! random number
           s1,         &
           s2,         &
           s12,        &
           z,          &
           tdm,        &   ! mean daily mean temperature
           trngm,      &   ! mean daily mean temperature
           tmaxm,      &   ! mean maximum temperature
           tminm,      &   ! mean minimum temperature
           tmaxd,      &   ! maximum temperatures for dry days
           tmaxw,      &   ! maximum temperatures for wet days
           tmaxe,      &   !'expected' maximum temperature for today
           tmine,      &   !'expected' minimum temperature for today
           tmaxs,      &   ! standard deviation in maximum temperature (K)
           tmins,      &   ! standard deviation in minimum temperature (K)
           cloudm,     &   ! mean cloud cover for today (fraction)
           cloudd,     &   ! dry day cloud cover
           cloudw,     &   ! wet day cloud cover
           cloude,     &   ! expected cloud cover today
           clouds,     &   ! standard deviation of cloud fraction
           v,          &
           tdum,       &   ! storage variable
           qdm,        &   ! mean relative humidity
           qdd,        &   ! dry day relative humidity
           qdw,        &   ! wet day relative humidity
           qde,        &   ! expected relative humidity (based on wet/dry decision)
           qdup,       &   ! upper bound of humidity distribution function
           qdlow,      &   ! lower bound of humidity distribution function
           y,          &
           b3,         &
           b2,         &
           b1,         &
           x1,         &
           amn,        &
           eud             ! expected daily average wind speed from monthly mean
!
      integer iyranom,    &
              iyrlast,    &
              nrun
!
      real precipfac,  &
           dif,        &
           humidfrac

#define WGEN_COMSAT
#include "inland_comsat.h"
!
! initialize this year's values of gdd
!
      if (iday.eq.1 .and. imonth.eq.1) then

!
!        call const (tcthis, npoi,  100.0)
        tcthis(:)=100.0
        twthis = -100.0
!
        gdd0this = 0.0
        gdd5this = 0.0
        gdd0cthis = 0.0
        gdd8this = 0.0
        gdd10this = 0.0
        gdd11this = 0.0
        gdd12this = 0.0
!
! initialize variables to zero at beginning of year
! for crop types that do not grow over two calendar years
!
! constant for gdd based on 10 cm soil temperature (planting bed-
! used to calculate leaf emergence after planting
!
          do  i = lbeg, lend
             do  j = 1, npft
!
              if (j .le. scpft-1) then  ! natural vegetation
                 ayanpp(i,j) = 0.0
              else if ( j .ge. scpft) then
                 ayanpp(i,j) = 0.0
              endif
	     enddo
         enddo
!
      endif
!
! initialize this crop year's values


    do  i = lbeg, lend
      do  j = scpft, ecpft

      if (iday.eq.pcd(j).and.imonth.eq.pcm(j)) then
!
	if(exist(i,13).eq.1.and.j.eq.13) then
           gddsoy(i,iyear)  = gddfzsoy(i)
           consdays(i)=0
           iniday(i)=9999
           maxcons(i)=0
           gsdays(i)=0
           gddfzcorn(i)=0.0
           gddfzsoy(i)=0.0
           gddfzsgc(i)=0.0
	else
	   gddsoy(i,iyear)=0.0
	endif

	if(exist(i,14).eq.1.and.j.eq.14) then
          gddcorn(i,iyear) = gddfzcorn(i)
          consdays(i)=0
          iniday(i)=9999
          maxcons(i)=0
          gsdays(i)=0
          gddfzcorn(i)=0.0
          gddfzsoy(i)=0.0
          gddfzsgc(i)=0.0
	else
          gddcorn(i,iyear)=0.0
	endif

	if(exist(i,16).eq.1.and.j.eq.16) then
          gddsgc(i,iyear)  = gddfzsgc(i)
          consdays(i)=0
          iniday(i)=9999
          maxcons(i)=0
          gsdays(i)=0
          gddfzcorn(i)=0.0
          gddfzsoy(i)=0.0
          gddfzsgc(i)=0.0
	else
	  gddsgc(i,iyear)=0.0
	endif

        if(croplive(i,j) .eq. 0 ) then
          gddplant(i,j)   = 0.0
          gddtsoi(i,j)    = 0.0
        endif

      endif

      enddo

    enddo
!
!
! ----------------------------------------------------------------------
! * * * set daily climatic variables for entire domain * * *
! ----------------------------------------------------------------------
!
      do 200 i = lbeg, lend
!
! ----------------------------------------------------------------------
! * * * other daily climate calculations * * *
! ----------------------------------------------------------------------
!
! calculated temperature extremes -- for vegetation limits (deg c)
!
! for this purpose, use the 10-day running mean temperature
!
        tcthis(i) = min (tcthis(i), (a10td(i) - 273.16))
        twthis(i) = max (twthis(i), (a10td(i) - 273.16))
!
! update this year's growing degree days
!

        gdd0this(i)  = gdd0this(i)  + max (0.0 ,(td(i) - 273.16))
        gdd0cthis(i) = gdd0cthis(i) + max (0.0 ,(td(i) - baset(15)))    ! wheat
        gdd5this(i)  = gdd5this(i)  + max (0.0 ,(td(i) - 278.16))
        gdd8this(i)  = gdd8this(i)  + max (0.0 ,(td(i) - baset(14)))    ! maize
        gdd10this(i) = gdd10this(i) + max (0.0 ,(td(i) - baset(13)))    ! soybean
        gdd11this(i) = gdd11this(i) + max (0.0 ,(td(i) - baset(17)))    ! oil palm
        gdd12this(i) = gdd12this(i) + max (0.0 ,(td(i) - baset(16)))    ! sugarcane


! form calculations of number of growing degree days between frost events
!
! events (e.g., tmin(i) .le. -2.2 C) this is applicable to CRU data
! differences exist when using a combination of CRU/NCEP
! -2.2 used as a threshold from Schwartz and Reiter, 2000.  International
! Journal of Climatology, 20: 929-932.  Changes in North American Spring
!

        if (tmin(i) .ge. 273.16) then
          consdays(i) = consdays(i) + 1
          maxcons(i)  = max(maxcons(i), consdays(i))

          if(maxcons(i) .eq. consdays(i)) then
            iniday(i) = cdays(i)+1 - maxcons(i)
          endif

           daygddc(1,cdays(i)) = min(30.0,max(0.0, (td(i) - baset(14))))
           daygdds(1,cdays(i)) = min(30.0,max(0.0, (td(i) - baset(13))))
           daygddsgc(1,cdays(i)) = min(30.0,max(0.0, (td(i) - baset(16))))
!
        else
           daygddc(1,cdays(i)) = 0
           daygdds(1,cdays(i)) = 0
           daygddsgc(1,cdays(i)) = 0
           consdays(i) = 0
        endif
!
         if (cdays(i).eq.365) then
!
! calculate total number of growing season GDD and number
! of days long
!
             if (iniday(i) .eq. 9999) then
                 iniday(i) = cdays (i)
                 maxcons(i) = 1
	     elseif (iniday(i) .eq. 0) then
	     iniday(i)=1
             endif
             endday(i) = iniday(i) + maxcons(i) - 1

          do 125 k = iniday(i), endday(i)

             gddfzcorn(i) =  gddfzcorn(i) + daygddc(1,k)
             gddfzsoy(i)  =  gddfzsoy(i)  + daygdds(1,k)

             gddfzsgc(i) =  gddfzsgc(i)   + daygddsgc(1,k)
             gsdays(i) = gsdays(i) + 1

125     continue

        elseif (cdays(i).eq.366) then

             if (iniday(i).eq.0) then
                 iniday(i) = 1
             endif

             endday(i) = iniday(i) + maxcons(i)-1

             gddfzcorn(i) = 0
             gddfzsoy(i)  = 0
             gddfzsgc(i) =  0
             gsdays(i) =    0

          do  k = iniday(i), endday(i)

             gddfzcorn(i) =  gddfzcorn(i) + daygddc(1,k)
             gddfzsoy(i)  =  gddfzsoy(i)  + daygdds(1,k)
             gddfzsgc(i) =  gddfzsgc(i)   + daygddsgc(1,k)
             gsdays(i) = gsdays(i) + 1
         enddo
        endif
!***************************** fim da provavel limpesa ****

!
! accumulate growing degree days for planted crops past planting
!
        do 150 j = scpft, ecpft
!
! for crops except winter wheat
! be careful with rotations
!
          if (croplive(i,j) .eq. 1.0 .and. iwheattype .ne. 2) then
!
                 gddplant(i,j) = gddplant(i,j) + max(0.0, min(td(i) &
                                 - baset(j), mxtmp(j)))
                 gddtsoi(i,j)  = gddtsoi(i,j) + max(0.0, min(tsoi(i,1) &
                                 - baset(j), mxtmp(j)))
!
! if winter wheat is planted, reduce thermal time accumulation
! by vernalization factor (calculated in crops.f) until crop
! is fully vernalized (cold-hardened)
!
             elseif (croplive(i,j) .eq. 1.0 .and. j .eq. 15 .and. &
                 iwheattype .eq. 2) then

                 gddplant(i,j) = gddplant(i,j) + vf(i) * max(0.0, min(td(i) &
                                 - baset(j), mxtmp(j)))
                 gddtsoi(i,j)  = gddtsoi(i,j)  +  vf(i) * max(0.0, min(tsoi(i,1) &
                                 - baset(j), mxtmp(j)))
          endif

	  if(j.eq.16) then

	     ik(i)=ik(i)+1
             if(cropy(i).eq.1.and.idpp(i,j).eq.1) then
	        print*,'GDD_start',i,iyear,jday
                gddpl15(i)=0
	        ik(i)=1
             endif

             if (ik(i).le.(mxmat(16)-30)) then
                gddpl15(i)= gddpl15(i) + max(0.0, min(td(i)- baset(j), mxtmp(j)))
             endif

             if (ik(i).eq.(mxmat(16)-30)) then
	        print*,'GDD_1/5 anos',i,iyear,jday,gddsgcp(i,1),gddpl15(i),mxmat(16)-30
                gddsgcp(i,1)= (gddsgcp(i,1)+ gddpl15(i))/2
	        print*,'GDD_1/5 anos',i,iyear,jday,gddsgcp(i,1)
	     endif

          endif

!
150  continue
!

200  continue
!
! return to main program
!
      return
      end subroutine dailymet
!
