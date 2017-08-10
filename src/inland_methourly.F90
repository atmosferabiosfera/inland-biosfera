#include "inland_config.h"
! ---------------------------------------------------------------------
      subroutine methourly (iyear,jday,timed,seed, seed2,seed3,seed4)
! ---------------------------------------------------------------------
!
! initializes surface meteorology from hourly data  
! avoids the need to go through much of subroutine 'daily' in weather.f 
!
!
      use inland_comwork
      use inland_comhyd
      use inland_comsum
      use inland_comhour
      use inland_comnitr
      use inland_parameters
      use inland_combcs
      use inland_comatm
      use inland_comcrop
      use inland_comsatp
!
      integer nhours,            &
              npts,              &
              seed,              & 
              seed2,             &
              seed3(32),             &
              seed4,             &
              jday,              &
              iyear,             &
              imonth,            &
              iday,              &
              ihtime,            &
              id1,               &
              id2,               &
              id3,               &
              n,                 &
              i,j,               &
              ihprecip,          &
              iyr,yfile,dfile,   &
              iyrinc
!      
      real*8 nr,tfile,        &   !nr - means that we are not reading these var
           rwork,                &
           var3sum,              &
           var3max,              &
           var3min,prec_30m_1h,  &
           chtime,               &
           psurfi,               &
           tdave,                &
           shum,                 &
           elevm, nada,          &
           ran2

      integer timed
!
      character*100 header
      character*100 metfile
      character*30  hstg
      character*39  pathn   ! has to be exact length of input string
      
!
      save nhours,          &
           var3sum,         &
           var3max,         &
           var3min
!
#define WGEN_COMSAT
#include "inland_comsat.h"

!         metfile = 'input/2008_2010.txt'
!csant- var1 to 6 are: prec, irradiance, air temp, relat. Hum., wind veloc., atm pressure
            n=1   !every time step reads the vars

 	    if(iyear.eq.imetyear.and.jday.eq.dmetyear.and.timed.eq.0) then

               open(232,file='input/single_agro.txt') !first time - open the file and let opened
! 2009      18       1   1.780  22.306   0.800   0.000 411.300  94.000
! ano      dia      hora vento  temp    prec    radin  rad-lonin    UR
! csant- var1 to 6 are: 1-prec, 2-irradiance,3-air temp, 4-relat. Hum., 5-wind veloc., 6-atm pressure


               do j = 1, nl  !find the time in the file

                  read(232,*) yfile,dfile,tfile, var5(n),var3(n), &
                              var1(n),var2(n),var6(n),var4(n)

                  if(iyear.eq.yfile.and.jday.eq.dfile.and.timed.eq.tfile*3600) goto 1212   ! find the date...

      	       enddo

	          write(*,*)'didnt find the time in the met inpunt file!'
	          stop

            else
	        
                prec_30m_1h= 0.
	     
                do j =1,1  !find the time, depens on IBIS dtime


                    read(232,*) yfile,dfile,tfile, var5(n),var3(n),var1(n), &
                                var2(n),var6(n),var4(n)

	            prec_30m_1h=  var1(n)           !prec_30m_1h + var1(n)

                    if(iyear.eq.yfile.and.jday.eq.dfile.and.timed.eq.tfile*3600) goto 1212   ! find the date...

                enddo

                write(*,*)'WARNING: model and file dates dont match!!'
                write(*,*),iyear,yfile,jday,dfile,timed,tfile*dtime,var5(n)
                stop
            endif

1212  continue

       var1(n) = prec_30m_1h
       rwork = (grav / rair / 0.0065)
!
! set elevation for the arlington grid cells
! based on elevation reading (meters) at morrisonville
!
       elevm = elevin !552.
!

         if (timed .eq. 0.0) then
              var3sum = 0.0
              nhours  =  0

              call const (var3max, npoi, 200.0)
	      call const (var3min, npoi, 400.0)
         endif              
!
              nhours  = nhours + 1

! csant- var1 to 6 are: prec, irradiance, air temp, relat. Hum., wind veloc., atm pressure

          var3min = min(var3min, var3(n) + 273.16)

          var3max = max(var3max, var3(n) + 273.16)

          psurfi  = 101325.0 * ((var3(n) + 273.16) / &
                    ((var3(n) + 273.16) + 0.0065 *   &
                    elevm)) ** rwork

          call const (psurf, npoi, psurfi)         ! surface pressure

!       call const (psurf, npoi, var6(n)*100)    ! reading surface pressure in hPa

          call const (precip, npoi, var1(n)* dtime / 3600.) ! precipitation (mm)  
	  call const (ta, npoi, var3(n)+ 273.16 )  ! air temperature
 
          var3sum = var3sum + var3(n)  + 273.16	   ! daily air temperature

          call const (cloud, npoi, var2(n))        ! insolation
!
          if (var4(n) .lt. 1.) then
              write(*,*) 'rh too low! ', var4(n)
              write(*,*) 'month, day, hour ', imonth, iday, ihtime,n
              stop 33
          endif 

! convert relative humidity to specific humidity
!
          shum  = var4(n) / 100. * qsat(esat(var3(n) + 273.16), &
                  psurfi)

!
          call const (qd, npoi, shum)             ! specific humidity 
          call const (ua, npoi, var5(n))           ! wind velocity
!

          if (timed .eq. (24.*3600.) - dtime) then

              call const (tmin, npoi, var3min)         ! minimum temp
              call const (tmax, npoi, var3max)         ! maximum temp

              var3sum = var3sum / float(nhours)

              tdave = var3sum
              call const (td, npoi, tdave)

              call dailymet(imonth, jday)
          endif

!
! return to main program
!
      return
      end
