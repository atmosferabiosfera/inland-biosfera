#include "inland_config.h"
module inland_combcs
      implicit none
      public
      save
!
!
! ------
! combcs
! ------
! xintopo(npoi): topography (m)
! xinveg(npoi): fixed vegetation map
! deltat(npoi): absolute minimum temperature - temp on average of coldest 
!              month (degree C)
      real*8, dimension(:), allocatable :: xintopo, xinveg, deltat

! xint(npoi,12): climatological temp + anomaly (C)
      real*8, dimension(:,:), allocatable :: xint

!---------------------------------------------------------------------------------
! Castanho HP, 2013 including initial variables input maps heterogeneous parameters

! xintauwood(npoi,npft):   initial tauwood map data (years)
! xinawood(npoi,npft):     initial awood map data (fraction)
! xinaroot(npoi,npft):     initial aroot map data (fraction)
! xinvmax(npoi,npft):      initial vmax rubisco map data (umol/m2/s)
! xinspecla(npoi,npft):     initial SLA map data (m2/g)

      real*8, dimension(:,:), allocatable :: xintauwood, xinawood, xinaroot, &
                                             xinvmax, xinspecla,xinaleaf
!-------------------------------------------------------------------------------
! land mask: 0: water; 1:land
      integer, dimension(:,:), allocatable :: lmask

! weather generator-specific variables
!--------------------------------------------------------
! FIXME: if these variables are used only by inland_daily subroutine, why they
!       need to make part of a module?
! xinwet(npoi,12): climatological wet days + anomaly (days/month)
! xinprec(npoi,12): climatological precipition + anomaly (mm/day)
! xintrng(npoi,12): climatological temp range + anomaly (C)
! xincld(npoi,12): climatological cloudiness + anomaly (%)
! xinq(npoi,12): climatological relative humidity + anomaly (%)
! xinwind(npoi.12): climatological wind speed + anomaly (m/s)
      real*8, dimension(:,:), allocatable :: xinwet, xinprec, xintrng, xincld, &
                                             xinq, xinwind

! xinprecd(npoi): daily climatological precipitation + anomaly * daily anomaly
!                 units: mm/day + mm/day * fraction
! xintd(npoi): daily climatological temperature + anomaly + daily anomaly
!              units: degrees C + degrees C + degrees C
! xintrngd(npoi): daily climatological temp range + anomaly * daily anomaly
!                 units: degrees C + degrees C + degrees C
! xintcldd(npoi): daily climatological cloud fraction + anomaly + daily anomaly
!                 units: fraction + fraction + fraction
! xinqd(npoi): daily climatological relative humidity + anomaly * daily anomaly
!              units: % + % * fraction
! xinwindd(npoi): daily climatological wind speed + anomaly * daily anomaly
!                 units: m/s + m/s * fraction
! xintmaxd(npoi): daily climatological maximum temperature
!                 units: degrees C 
! xintmind(npoi): daily climatological minimum temperature
!                 units: degrees C 

      real*8, dimension(:), allocatable :: xinprecd,xintd,xintrngd,xincldd, &
                                           xinqd,xinwindd,xintmaxd,xintmind

      real*8, dimension(:,:), allocatable :: xinwetmon, xinprecmon, xintrngmon, xincldmon, &
                                             xinqmon, xinwindmon, xintmon

! Anomaly-specific variables
!-----------------------------
! clmt(npoi,12):    climatological temperature (C)
! clmprec(npoi,12): climatological precipitation (mm/day)
! clmcld(npoi,12):  climatological cloudiness (%)
! clmq(npoi,12):    climatological relative humidity (%)
      real*8, dimension(:,:), allocatable :: clmt,clmprec,clmcld,clmq

! obst(npoi,12):    observed temperature (C)
! obsprec(npoi,12): observed precipitation (mm/day)
! obscld(npoi,12):  observed cloudiness (%)
! obsq(npoi,12):    observed relative humidity (%)
! obswind(npoi,12):    observed wind (m/s)
      real*8, dimension(:,:), allocatable :: obst,obsprec,obscld,obsq,obswind

! obswet(npoi,12):  observed wet days (days/month)
! obstrng(npoi,12): observed temperature range (degrees C)
      real*8, dimension(:,:), allocatable :: obswet,obstrng
! Readit-specific variables (??)
! clmwet(npoi,12):  climatological wet days (days/month)
! clmtrng(npoi,12): climatological temperature range (degrees C)
      real*8, dimension(:,:), allocatable :: clmwet,clmtrng,clmwind
!
!              !Station Data         	
!       stintd(npoi),      ! daily Station temperature (C) 
!       stinqd(npoi),      ! daily Station relative humidity (%) 
!       stinprecd(npoi),   ! daily Station precipitation (mm/day) 
!       stinwindd(npoi),   ! daily Station windspeed (m/s) 
!       stincldd(npoi),    ! daily Station cloud fraction (fraction) 
!       stintmax(npoi),    ! daily Station 2m temp maximum (C)  
!       stintmin(npoi),    ! daily Station 2m temp minimum (C)
     real*8, dimension(:), allocatable :: stintd, stinqd,  &
                                          stinprecd, stinwindd, stincldd,&
                                          stintmax, stintmin
! FIXME: We have some confusion here. On some places, 'mpt' is used as
!       alternative to npoi; other places have 'lbeg:lend' for the same reason.
end module inland_combcs
