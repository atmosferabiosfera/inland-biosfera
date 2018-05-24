#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine drystress(kpti, kptj)
! ---------------------------------------------------------------------
      use inland_parameters, only: nsoilay, stressfac, dtime, npoi, epsilon
      use inland_control
      use inland_comveg
      use inland_comsoi
      use inland_comcrop, only:isimagro
      use inland_comhour, only:imetyear

      implicit none

!-----------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer i, k, kk

      real*8 awc,   		& ! available water content (fraction)
             znorm, 		& ! normalizing factor
             zwilt,         & ! function of awc, =1 if awc = 1 (no stress)
             lamdau,   	     & ! upper layer parameter determining the relative importance of the roots
             lamdal,   	     & ! low layer
             uprate(nsoilay),& ! maximum water uptake rate
             aw(nsoilay),    & ! awailable water fraction for six layer
             awc1,      	 & ! awailable water fraction for 1/3 upper rooting zone 
             awc2,      	 & ! awailable water fraction for the depth from 1/3 root zone to bottom of root zone 
             updep,     	 & ! depth of the upper 1/3 root zone
             lowdep,    	 & ! depth of lower 2/3 root zone
             depf1, depf2      ! fraction

! stressfac determines the 'strength' of the soil moisture
! stress on physiological processes
!
! strictly speaking, stresst* is multiplied to the vmax
! parameters used in the photosynthesis calculations
!
! stressfac here determines the shape of the soil moisture response
      znorm = 1.0 - exp(stressfac)
      do 100 i = kpti, kptj
! initialize stress parameter
         stresstl(i) = 0.0
         stresstu(i) = 0.0

!isimrwu = 0  , ! isimrwu - Root water uptake module - 0: according to Foley et al., 1996; 1 according to Li et al. (2006) (default 0)

         if (isimrwu.eq.1) then    ! Castanho Kai RootWaterUptake
! 
! Kaiyuan Li

           stre_tl(i) = 0.0
           stre_tu(i) = 0.0

! ----------------------------------------------------------------
! added by Kaiyuan Li to calculate lamda values
!
           do k = 1, nsoilay 
              aw(k) = min (1.0, max (dble(0.0),                      &
                      (wsoi(i,k)*(1 - wisoi(i,k)) - swilt(i,k)) /   &
                      (sfield(i,k) - swilt(i,k)) ) )
           end do
!      
! initialize awc1, awc2, updep, lowdep
          awc1 = 0
          awc2 = 0
          updep = 0
          lowdep = 0


! Lower canopy: calculating available water content for upper and lower depth  
          do k = 1, nsoilay
             updep = updep + 100*hsoi(k)
             if (abs(updep - drl/3.0) .le. 2.0) then
                awc1 = awc1 + aw(k)*hsoi(k)*100/(drl/3.0)
                do kk = k+1, nsoilay
                   lowdep = lowdep + 100*hsoi(kk)
                   if(abs(lowdep - 2.0*drl/3.0) .le. 2.0) then
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*drl/3.0)
                      go to 88
                   else if ((lowdep - 2.0*drl/3.0) .gt. 2.0) then
                      depf2 = (lowdep - 2.0*drl/3.0)/(hsoi(kk)*100)
                      awc2 = awc2 + aw(kk)*(1.0 - depf2)*hsoi(kk)*100/(2.0*drl/3.0)
                      go to 88
                   else
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*drl/3.0)
                   end if
                end do
             else if((updep - drl/3.0) .gt. 2.0) then
                depf1 = (updep - drl/3.0)/(hsoi(k)*100)
                awc1 = awc1 + aw(k)*(1.0 - depf1)*hsoi(k)*100/(drl/3.0)
                do kk = k+1, nsoilay
                   lowdep = lowdep + 100*hsoi(kk)
		           if(abs(lowdep - 2.0*drl/3.0) .le. 2.0) then
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*drl/3.0)
                      go to 88
                   else if ((lowdep - 2.0*drl/3.0) .gt. 2.0) then
                      depf2 = (lowdep - 2.0*drl/3.0)/(hsoi(kk)*100)
                      awc2 = awc2 + aw(kk)*(1.0 - depf2)*hsoi(kk)*100/(2.0*drl/3.0)
                      awc2 = awc2 + depf1*hsoi(k)*100/(2*drl/3.0)
                      go to 88
                   else
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*drl/3.0)
                   end if                
                end do
             else
                awc1 = awc1 + aw(k)*hsoi(k)*100/(drl/3.0)
             endif
          end do       
  
88        if((awc1 .ge. 0.2 .and. awc1 .le. 0.5) .and. (awc2 .gt. 0.5)) then 
             lamdal = 0.75
          else if ((awc1 .lt. 0.2) .and. (awc2 .gt. 0.5)) then 
             lamdal = 0.50
          else if ((awc1 .lt. 0.2) .and. (awc2 .ge. 0.2 .and. awc2 .le. 0.5)) then
             lamdal = 0.75
          else if ((awc1 .gt. 0.5) .and. (awc2 .ge. 0.2 .and. awc2 .le. 0.5)) then
             lamdal = 1.25
          else if ((awc1 .gt. 0.5) .and. (awc2 .lt. 0.2)) then
             lamdal = 1.5
          else if ((awc1 .ge. 0.2 .and. awc1 .le. 0.5) .and. (awc2 .lt. 0.2)) then
             lamdal = 1.25
          else
             lamdal = 1.00  
          end if

! initialize awc1, awc2, updep, lowdep
          awc1 = 0
          awc2 = 0
          updep = 0
          lowdep = 0

! Upper canopy: calculating available water content for upper and lower depth
          do k = 1, nsoilay
             updep = updep + 100*hsoi(k)
             if (abs(updep - dru/3.0) .le. 2.0) then
                awc1 = awc1 + aw(k)*hsoi(k)*100/(dru/3.0)
                do kk = k+1, nsoilay
                   lowdep = lowdep + 100*hsoi(kk)
                   if(abs(lowdep - 2.0*dru/3.0) .le. 2.0) then
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*dru/3.0)
                      go to 99
                   else if ((lowdep - 2.0*dru/3.0) .gt. 2.0) then
                      depf2 = (lowdep - 2.0*dru/3.0)/(hsoi(kk)*100)
                      awc2 = awc2 + aw(kk)*(1.0 - depf2)*hsoi(kk)*100/(2.0*dru/3.0)
                      go to 99
                   else
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*dru/3.0)
                   end if
                end do
             else if((updep - dru/3.0) .gt. 2.0) then
                depf1 = (updep - dru/3.0)/(hsoi(k)*100)
                awc1 = awc1 + aw(k)*(1.0 - depf1)*hsoi(k)*100/(dru/3.0)
                do kk = k+1, nsoilay
                   lowdep = lowdep + 100*hsoi(kk)
                   if(abs(lowdep - 2.0*dru/3.0) .le. 2.0) then
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*dru/3.0)
                      go to 99
                   else if ((lowdep - 2.0*dru/3.0) .gt. 2.0) then
                      depf2 = (lowdep - 2.0*dru/3.0)/(hsoi(kk)*100)
                      awc2 = awc2 + aw(kk)*(1.0 - depf2)*hsoi(kk)*100/(2.0*dru/3.0)
                      awc2 = awc2 + depf1*hsoi(k)*100/(2*dru/3.0)
                      go to 99
                   else
                      awc2 = awc2 + aw(kk)*hsoi(kk)*100/(2.0*dru/3.0)
                   end if
                end do
             else
                awc1 = awc1 + aw(k)*hsoi(k)*100/(dru/3.0)
             endif
          end do

99        if((awc1 .ge. 0.2 .and. awc1 .le. 0.5) .and. (awc2 .gt. 0.5)) then
             lamdau = 0.75
          else if ((awc1 .lt. 0.2) .and. (awc2 .gt. 0.5)) then
             lamdau = 0.50
          else if ((awc1 .lt. 0.2) .and. (awc2 .ge. 0.2 .and. awc2 .le. 0.5)) then
             lamdau = 0.75
          else if ((awc1 .gt. 0.5) .and. (awc2 .ge. 0.2 .and. awc2 .le. 0.5)) then
             lamdau = 1.25
          else if ((awc1 .gt. 0.5) .and. (awc2 .lt. 0.2)) then
             lamdau = 1.5
          else if ((awc1 .ge. 0.2 .and. awc1 .le. 0.5) .and. (awc2 .lt. 0.2)) then
             lamdau = 1.25
          else
             lamdau = 1.00
          end if
! --------------------------------------------------------------
! ----- Kaiyuan Li to change to origianl root water uptake model
!       lamdau = 1.0
!       lamdal = 1.0
! --------------------------------------------------------------
!
! fraction of soil water uptake in each layer
!
        do 110 k = 1, nsoilay
!
! plant available water content (fraction)
!
          awc = min (1.0, max (dble(0.0),                 &
               (wsoi(i,k)*(1 - wisoi(i,k)) - swilt(i,k)) /&
               (sfield(i,k) - swilt(i,k))    ))

          zwilt = (1. - exp(stressfac * awc)) / znorm

! update for each layer
!
! ------------------------------------------------------------------
! added by Kaiyuan Li to calculate Maximum Water Uptake Rate (uprate)
! see Campbell and Norman 1998, An introduction to environmental biophysics
          uprate(k) = 1.0 - (1.0 + 1.3 * awc)**(-5.0)  
! -------------------------------------------------------------------
! modified by Kaiyuan Li for water-uptake 
!
          stressl(i,k) = (froot(k,1)**lamdal) * max (dble(0.0), min (1.0, zwilt))   ! added by Kaiyuan Li
          stressu(i,k) = (froot(k,2)**lamdau) * max (dble(0.0), min (1.0, zwilt))   ! added by Kaiyuan Li
!
! integral over rooting profile
!
           stre_tl(i) = stre_tl(i) + stressl(i, k)  !added by Kaiyuan Li
           stre_tu(i) = stre_tu(i) + stressu(i, k)  !added by Kaiyuan Li 

110    continue

        do k = 1, nsoilay  ! added by Kaiyuan Li
          stresstl(i) = stresstl(i) + (stressl(i,k)/max(stre_tl(i), epsilon)) * &
                       max(0.0, min(1.0, uprate(k))) !  by Kaiyuan Li
          stresstu(i) = stresstu(i) + (stressu(i,k)/max(stre_tu(i), epsilon)) * &
                       max(0.0, min(1.0, uprate(k)))  ! by Kaiyuan Li
        end do ! added by Kaiyuan Li

! Castanho Kai end--------------------------------------------------------------------------------------------------
         
             else          ! Castanho Kai 

! fraction of soil water uptake in each layer
         do 111 k = 1, nsoilay

! plant available water content (fraction)
            awc = min(dble(1.0),max(dble(0.0),(wsoi(i,k)*(1 - wisoi(i,k)) - swilt(i,k)) / &
                      (sfield(i,k) - swilt(i,k))))

         if(isimagro .eq. 0) then
            zwilt = (1. - exp(stressfac * awc)) / znorm
         else
            zwilt = dble(1.0) - (log(1+799.0*exp(-38.0*awc))/alog(800.))
         endif

! update for each layer
            stressl(i,k) = froot(k,1) * max (dble(0.0), min (dble(1.0), zwilt))
            stressu(i,k) = froot(k,2) * max (dble(0.0), min (dble(1.0), zwilt))

! integral over rooting profile
            stresstl(i) = stresstl(i) + stressl(i,k)
            stresstu(i) = stresstu(i) + stressu(i,k)
111      continue

            endif          ! Castanho Kai RootWaterUptake
100   continue

!      if(dtime.eq.39600.and.(iyear.eq.2010.and.jday.le.120))then
      if(imetyear .ne. 9999)then
            write(225,34)jday,wsoi(1,1),wsoi(1,2),wsoi(1,3),wsoi(1,4),wsoi(1,5),wsoi(1,6),  &
                         wsoi(1,7),wsoi(1,8),wsoi(1,9),wsoi(1,10),wsoi(1,11),               &
                         sfield(1,1),swilt(1,1),stresstl(1)

 34   format (1x,i3,1x,15(1x,f4.2))
      
      endif 


      return
end subroutine drystress
