#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine irrigation(ilens, irrigate)
! ---------------------------------------------------------------------
!
! routine that calculates the amount of water that is applied to a
! managed ecosystem on a daily basis (mm/day)
!
! based on average daily water content in the soil
!
! this amount will be evenly applied in timestep increments (start and stop time)
! based on the duration of the event - similar to the way precipitation
! events are handled 
!
! common blocks
!
!
      use inland_parameters
      use inland_control, only:iday, imonth
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
!
! local variables
!
      integer i, j, k
!
      real*8 awcmax,  &
             awc 

!
      do 100 i = lbeg, lend
         if (iday .eq. 1 .and. imonth .eq. 1) then
            totirrig(i) = 0.0
         endif
         awcmax = 0.0
         awc    = 0.0 
!
         do 110 k = 1, 5  ! top 5 soil layers = 100 cm 
!
! calculate the maximum available water in top 1 m (root zone) at field capacity
!
            awcmax = awcmax + max(0.0, (sfield(i,k) - swilt(i,k))) * &
                             hsoi(k) * poros(i,k) * 100
!
! calculate actual amount of water available to plant - current water content
! based on daily average water and ice content in soil layers down to a meter 
! 
            awc    = awc    + max(0.0, (adwisoilay(i,k)   + &
                            (1. - adwisoilay(i,k))        * &
                            adwsoilay(i,k)) - swilt(i,k)) * &
                            hsoi(k) * poros(i,k)          * 100
!
110      continue    
!
!    irrigation will occur if :
!
!  * the crop has been planted
!  * the minimum daily temperature is > 5 C
!  * the 5-day running mean temperature is > 10 C 
!  * the actual soil water content in the top 1 m of soil
!    is less than or equal to 50% of the maximum value at field capacity,
!
         if (awc .le. 0.50 * awcmax &
           .and. (tmin(i) .gt. 278.16 .and. a5td(i) .gt. 283.16) &
           .and. (croptype(i) .ge. scpft)) then
!           .and. (croplive(i,j) .gt. 0)  ! have to add in that irrigation only occurs when crop is planted! 
!
! irrigation (applied above the canopy)  is used to make up the difference to field capacity
! convert awc values from cm to mm - this is a per day value 
! so it is consistent with the values used in weather.f (diurnal) for precipitation 
!
! set upper limit based on literature search typical of what a farmer could
! apply in a typical day of irrigation and what rates of application are
! 
            xirrig(i) = min(150.0, max(0.0, awcmax - awc) * 10.0)
         else
            xirrig(i) = 0.0
        endif
!
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
        if (time.ge.starti .and. time.lt.endi & 
           .and. precip(i) .eq. 0.00) then  
!
         xirriga(i) = xirrig(i) / ilens

!
! update annual total - totirrig
! rate of irrigation multiplied by length of timestep (mm/s * s) = total applied
! for this timestep 
!
!          totirrig(i) = totirrig(i) + (xirriga(i) * dtime)
!
        endif

100   continue
!
! return to main program
!
      return
!      
end subroutine irrigation
