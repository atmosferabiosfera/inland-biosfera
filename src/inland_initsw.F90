#include "inland_config.h"
! -----------------------------------------------------------------------------------      
subroutine initsw(kpti,kptj)         
! -----------------------------------------------------------------------------------
! Castanho Kaiyuan Li for Green-Ampt infiltration  NOTE NEW FUNCTION
!
! ---- To save soil moisture content before precipitation event so that
!      the Greeb-Ampt is able to calculate the change in volumetric soil water content
!      in the wetting depth
! common blocks
       use inland_parameters
       use inland_comatm
       use inland_comsoi

       implicit none

! input variables
      integer kpti        ! index of 1st point of little vector in big lpt vector
      integer kptj,i,j    ! index of last point of little vector Castanho Kai NOTE included ie j not included???
!       
      do i = kpti, kptj  !Castanho Kai NOTE changed to kpti kptj describe variables!!
          if((iwet(i) .eq. 1) .and. (aint(t_sec/dtime) .eq. aint(t_startp/dtime)) &
            .and. (wpud(i) .lt. 0.000001)) then
             do j = 1, nsoilay
                swater(i, j) = wsoi(i, j)
                sice(i, j) = wisoi(i, j)
             end do
          end if
       end do
!
end subroutine initsw         
