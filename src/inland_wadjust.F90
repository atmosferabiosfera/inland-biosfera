#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine wadjust(kpti, kptj)
! ---------------------------------------------------------------------
! set wsoi, wisoi to exactly 0 if differ by negligible amount, 
! to protect epsilon logic in soilctl and soilh2o
!
! ice-liquid transformations in soilctl loop 400 can produce very
! small -ve amounts due to roundoff error, and very small -ve or +ve
! amounts can cause (harmless) "underflow" fpes in soilh2o
! ---------------------------------------------------------------------
      use inland_comsoi
      use inland_com1d
      use inland_parameters
      use inland_comhyd

      implicit none
!-----------------------------------------------------------------------
! input variable
      integer kpti            ! index of 1st point of little vector
                             ! in big lpt vector
      integer kptj            ! index of last point of little vector
 
! local variables
      integer k, i
      real*8 ztot0, ztot1
!-----------------------------------------------------------------------

      do 100 k = 1, nsoilay
         do 110 i = kpti, kptj

! initial total soil water
            ztot0 = hsoi(k) * poros(i,k) * rhow * ((1. - wisoi(i,k)) * &
                    wsoi(i,k) + wisoi(i,k))

! set bounds on wsoi and wisoi
            if (wsoi(i,k).lt.epsilon) wsoi(i,k)  = 0.0
            if (wisoi(i,k).lt.epsilon) wisoi(i,k) = 0.0
            wsoi(i,k)  = min (dble(1.), wsoi(i,k))
            wisoi(i,k) = min (dble(1.), wisoi(i,k))
            if (wisoi(i,k).ge.1-epsilon) wsoi(i,k) = 0.0

! for diagnosis of total adjustment
            ztot1 = hsoi(k) * poros(i,k) * rhow * ((1. - wisoi(i,k)) * &
                    wsoi(i,k) + wisoi(i,k))
            gadjust(i) = gadjust(i) + (ztot1 - ztot0) / dtime
110      continue
100   continue
      return
end subroutine wadjust
