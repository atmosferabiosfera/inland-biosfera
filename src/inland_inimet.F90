#include "inland_config.h"
! ---------------------------------------------------------------------
      subroutine inimet (iyear, imonth, iday, nday, time, seed)
! ---------------------------------------------------------------------
!
!   initializes surface meteorology from hourly data  
! avoids the need to go through much of subroutine 'daily' in weather.f 
!
!
      use inland_parameters
      use inland_comwork
      use inland_comatm
      use inland_comhyd
      use inland_comsum
      use inland_combcs
      use inland_comhour
      use inland_comcrop
      use inland_comnitr
      use inland_comsatp
!
!
	integer nhours,npts,seed,nday,iyear,imonth,iday,&
		ihtime,id1,id2,id3,n,i,ihprecip,iyr,iyrinc
!
!      
	real*8 time,rwork,var3sum,var3max,var3min,chtime,psurfi,tdave,&
	       shum,elevm,ran2 
!
	character*100 header
	character*100 metfile
	character*30  hstg
	character*39  pathn   ! has to be exact length of input string
      
!
      save nhours,var3sum,var3max,var3min
!
#define WGEN_COMSAT
#include "inland_comsat.h"  
!
      pathn = '/Volumes/wet/kucharik/IBIS_input/input/'
!
!   spin up period by choosing a random year
! of hourly weather from 1986-2000
! have to make sure that npts matches the number
! of hourly points in a particular year (iyear)
! during spin up
!
      npts = 8760  ! number of yourly points in non-leap year 
!
!   generate the hourly met file to use (iyr)
! spin-up
!
      if (iyear .lt. 1986 .or. iyear .gt. 2000) then
!
!   is the actual calendar year (iyear) a leap year?
! if so, change the number of points to read in
!
        if (mod(iyear,4) .eq. 0) then
          if (mod(iyear,100) .ne.0) then
            npts = 8784
          else if (mod(iyear/100,4).eq.0) then
            npts = 8784
          endif
        endif
!
!   then choose a leap year out of met record  
! but only from 92, 96, or 2000 since 88 was
! an anomoly year (drought)
!
        if (npts .eq. 8784) then 
           iyrinc = min(8, int(ran2(seed) * 9))
           iyrinc = int(iyrinc/4) * 4
           iyr    = 1992 + iyrinc 
!
!   not a calendar leap year - choose any year
! choose any other year - but note, if you were
! to choose a leap year, it will only read in the
! first 8760 pts from the met file
!
        else
          iyrinc = min(14, int(ran2(seed) * 15)) 
          iyrinc = max(0, iyrinc)
          iyr    = 1986 + iyrinc 
!
        endif
!
      else
!
! now met files match calendar year 
!
! change npts for leap years
!
        if (iyear .eq. 1988 .or. &
           iyear .eq. 1992 .or. &
           iyear .eq. 1996 .or. &
           iyear .eq. 2000) npts = 8784             
!
        iyr = iyear
!
      endif
!

! arlington met station data files
!------------------------------------------
      if (iyr .eq. 1986) then 
         metfile = pathn//'arlington/metdata/1986'
      endif
!
      if (iyr .eq. 1987) then
         metfile = pathn//'arlington/metdata/1987'
      endif
!
      if (iyr .eq. 1988) then
         metfile = pathn//'arlington/metdata/1988'
      endif
!
      if (iyr .eq. 1989) then 
         metfile = pathn//'arlington/metdata/1989'
      endif
!
      if (iyr .eq. 1990) then 
         metfile = pathn//'arlington/metdata/1990'
      endif
!
      if (iyr .eq. 1991) then 
         metfile = pathn//'arlington/metdata/1991'
      endif
!
      if (iyr .eq. 1992) then
         metfile = pathn//'arlington/metdata/1992_lancaster'
      endif
!
      if (iyr .eq. 1993) then
         metfile = pathn//'arlington/metdata/1993'
      endif
!
      if (iyr .eq. 1994) then
         metfile = pathn//'arlington/metdata/1994'
      endif
!
      if (iyr .eq. 1995) then
         metfile = pathn//'arlington/metdata/1995'
      endif
!
      if (iyr .eq. 1996) then
         metfile = pathn//'arlington/metdata/1996'
      endif
!
      if (iyr .eq. 1997) then
         metfile = pathn//'arlington/metdata/1997'
      endif
!
      if (iyr .eq. 1998) then
         metfile = pathn//'arlington/metdata/1998'
      endif
!
      if (iyr .eq. 1999) then
         metfile = pathn//'arlington/metdata/1999'
      endif
!
      if (iyr .eq. 2000) then
         metfile = pathn//'arlington/metdata/2000'
      endif
!
!-----------------------------------------------
      if (imonth .eq. 1 .and. iday .eq. 1 .and. &
         time .eq. 0.0) then
        write(6,*) 'Reading hourly weather data ...'
        write(6,*) 'Using met file year = ', iyr
        open(12,file = metfile)
        read(12,*) header
        do 100 n = 1, npts 
          read(12,500) iiyear(n),iimonth(n),iiday(n), iihour(n)
 100    continue
 500    format(i4,1x,i2,1x,i2,1x,i4)
!
        close(12)    
!
        open(12,file = metfile)
        read(12,*) header
        do 200 n = 1, npts 
          read(12,*) hstg,id1,var1(n),var2(n),var3(n), &
                            var4(n),var5(n),var6(n)
 200    continue
        close(12)
!
      endif
!
        rwork = (grav / rair / 0.0065)
!
!   set elevation for the arlington grid cells
! based on elevation reading (meters) at morrisonville
!
        elevm = 293.8
!
        do 400 n  = 1, npts 
!
!         if (iyear .eq. iiyear(n)) then
            if (imonth .eq. iimonth(n)) then
              if (iday .eq. iiday(n))     then
!
                 chtime = (time / 3600.)
                 ihtime = int(chtime) + 1
!
                if (ihtime .eq. (iihour(n) / 100)) then
                  if (time .eq. 0.) then
                    var3sum = 0.0
                    nhours  = 0
                    call const (tmin, npoi, 400.0)
                    call const (tmax, npoi, 200.0)
                  endif              
!
                  nhours  = nhours + 1
                  var3min = min(tmin(1), var3(n) + 273.16)
                  var3max = max(tmax(1), var3(n) + 273.16)
                  psurfi  = 101325.0 * ((var3(n) + 273.16) / &
                           ((var3(n) + 273.16) + 0.0065 *    &
                             elevm)) ** rwork

                  call const (psurf, npoi, psurfi)         ! surface pressure
!
!   for precipitation, value at time n is precip accumulated since last
! time step, so we need to get the precip value for the time step n+1 
!
                  ihprecip = min(n+1, npts)  ! use precip from the next time step
                  call const (precip, npoi, var1(ihprecip) &
                            * dtime / 3600.)             ! precipitation (mm)  
                  call const (tmin, npoi, var3min)         ! minimum temp
                  call const (tmax, npoi, var3max)         ! maximum temp
                  call const (ta, npoi, var3(n) + 273.16)  ! air temperature
                  var3sum = var3sum + var3(n) + 273.16
                  call const (cloud, npoi, var2(n))        ! insolation
!
                  if (var4(n) .lt. 1.) then
                    print *, 'rh too low! ', var4(n)
                    print *, 'month, day, hour ', imonth, iday, ihtime,n
                    stop 33
                  endif 
!               
! convert relative humidity to specific humidity
!
                  
              shum  = var4(n) / 100. * qsat(esat(var3(n) + 273.16), &
                         psurfi) 
!
                  call const (qd, npoi, shum)             ! specific humidity 
                  call const (ua, npoi, var5(n))           ! wind velocity
!
                  if (time .eq. 24. * 3600. - dtime) then
                    if (nhours .ne. 24. * int (3600./dtime)) then
                      print *, 'this day has incorrect hours!'
                      print *, 'nhours = ', nhours
                      print *, 'iyear ', iyr 
                      stop 7
                    endif
                    var3sum = var3sum / float(nhours)
!
! three ways to compute daily average temp
! * average of hours        :  tdave = var3sum
! * weighted of min and max :  tdave = 0.44 * var3max + 0.56 * var3min
! * average of min and max  :  tdave = (var3min + var3max) / 2 
!
                   tdave = var3sum
                   call const (td, npoi, tdave)
!                  call dailymet(imonth, iday, seed, tdave, 0)
                   call dailymet(imonth, iday, seed, 1)
!
                 endif     
                endif
              endif
            endif
!         endif
!
 400    continue
!
!
! return to main program
!
      return
end subroutine inimet 
