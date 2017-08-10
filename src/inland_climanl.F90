#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
! this subsroutine is only used to initialize growing degree days,
! coldest temp, and warmest temp at very beginning - provides a
! climate 'history' based on monthly mean values
subroutine climanl
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comveg
      use inland_combcs
      use inland_comcrop

      implicit none

! local variables
      integer it1w,  &    ! indice of previous month (interpolation)
              it2w,  &    ! indice of following month (interpolation)
              i,k,lda     ! loop indices

      real*8 rwork, &     ! work variable (1/ndaypm)
             dt,    &     ! used for interpolation
             dtemp        ! interpolated temperature

! initialize values
      do 100 i = lbeg, lend

! coldest monthly temperature (year 0) in deg c
         tc(i) = min (xint(i,1),  xint(i,2),  xint(i,3),  &
                      xint(i,4),  xint(i,5),  xint(i,6),  &
                      xint(i,7),  xint(i,8),  xint(i,9),  &
                      xint(i,10), xint(i,11), xint(i,12))

! warmest monthly temperature (year 0) in deg c
         tw(i) = max (xint(i,1),  xint(i,2),  xint(i,3),  &
                      xint(i,4),  xint(i,5),  xint(i,6),  &
                      xint(i,7),  xint(i,8),  xint(i,9),  &
                      xint(i,10), xint(i,11), xint(i,12))

         tcmin(i) = tc(i) + deltat(i)
100   continue 

! interpolating climatological monthly input values to daily
      do 200 i = lbeg, lend 
         gdd0(i) = 0.
         gdd5(i) = 0.
         do 210 k = 1, 12
            rwork = 1. / float(ndaypm(k))
            do 220 lda = 1, ndaypm(k)
               if (float(lda).lt.float(ndaypm(k)+1)*0.5) then
                  it1w = k - 1
                  it2w = k
                  dt   = (float(lda) - 0.5) * rwork + 0.5
               else
                  it1w = k
                  it2w = k + 1
                  dt   = (float(lda) - 0.5) * rwork - 0.5
               end if

               if (it1w.lt. 1) it1w = 12
               if (it2w.gt.12) it2w = 1

               dtemp = xint(i,it1w) + dt * (xint(i,it2w) - xint(i,it1w))

! growing degree days, using deg c
               gdd0(i) = gdd0(i) + max(dble(0.0), dtemp)
               gdd5(i) = gdd5(i) + max(dble(0.0), (dtemp - 5.0))


     if(isimagro .gt. 0)then
        gdd0c(i) = gdd0c(i) + max(0.0, min(dtemp  +273.16-baset(15),  26.0))
        gdd8(i)  = gdd8(i)  + max(0.0, min(dtemp  +273.16-baset(14),  30.0))
        gdd10(i) = gdd10(i) + max(0.0, min(dtemp  +273.16-baset(13), 30.0))
        gdd12(i) = gdd12(i) + max(0.0000001, min(dtemp  +273.16-baset(16), 30.0))

	

	    if(pmmin(16).lt.pcm(16)-2) then

	       if(k.ge.pmmin(16).and.k.lt.pcm(16)-2) then  
           gddsgcp(i,1)= gddsgcp(i,1) + max(0.0, min(dtemp+273.16-baset(16),30.0))
           endif


	    else if(pmmin(16).gt.pcm(16)-2) then

	        if(k.ge.pmmin(16).and.k.le.12) then
               gddsgcp(i,1)= gddsgcp(i,1) + max(0.0, min(dtemp+273.16-baset(16),30.0))
            endif

            if(k.ge.1.and.k.lt.pcm(16)-2) then
              gddsgcp(i,1)= gddsgcp(i,1) + max(0.0, min(dtemp+273.16-baset(16),30.0))
	        endif

	    endif

       endif ! check for crop existence
   
 220      continue
 210    continue
     if(isimagro .gt. 0)then
	   gddsgcp(i,2) = gdd12(i)
	   gddsgcp(i,1)= gddsgcp(i,1) + gddsgcp(i,2)
     endif

 200  continue

! call routine to determine pft existence arrays
      call existence

      return
end subroutine climanl
