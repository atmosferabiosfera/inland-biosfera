#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine coldstart(rdlsf, lati, loni)
! ---------------------------------------------------------------------
!     
!     initialize some model varaibles for cold start conditions
!     
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsoi
      use inland_comsno
      use inland_combcs
      use inland_comcrop, only:isimagro

      implicit none
!-----------------------------Arguments---------------------------------
! Input arguments
!
      logical rdlsf           ! true: initial data for Antarctica, Greenland.
!
      real*8 lati(npoi)          ! latitude of land point on long vector [rads]
      real*8 loni(npoi)          ! longitude of land point on long vector [rads]
      integer k, l
!
   if(isimagro .gt. 0) then

      do k = lbeg, lend
         fi(k) = 0.
      enddo

      do l = 1, nsnolay
      do k =lbeg, lend
         hsno(k,l) = 0.
         tsno(k,l) = 273.16
      enddo
      enddo

      do l = 1, nsoilay
         do k =lbeg, lend
            tsoi(k,l) = 278.16
            wsoi(k,l) = 0.50
            wisoi(k,l) = 0.
         enddo
      enddo
   else

      fi(:) = 0.
      hsno(:,:) = 0.
      tsno(:,:) = 273.16
      wisoi(:,:) = 0.

#ifndef SINGLE_POINT_MODEL
      do l = lbeg, lend
         tsoi(l,:) = 293.16 + 15.*cos(2.*lati(l)*pi/180.)
      enddo
#else /* SINGLE_POINT_MODEL */
!      tsoi(:,:) = 278.16
#endif /* SINGLE_POINT_MODEL */
    endif ! ! check for crop existence

      do k =lbeg, lend
         tg(k) = tsoi(k,1)
      enddo

! Initialize temperature and snow depths in Antarctica and Greenland
      if (rdlsf) then 
         do k=lbeg, lend

! Antarctica
            if (lati(k) .le. -60.) then
               do l = 1, nsoilay
                  tsoi(k,l) = xint(k,1) + 273.16
                  wisoi(k,l) = 1.0
               end do
               tg(k) = tsoi(k,1)
               hsno(k,1) = 0.05
               hsno(k,2) = 2.
               hsno(k,3) = 2.
               fi(k) = 1.
            end if

! Greenland
            if (lati(k) .ge. 60  .and. &
                lati(k) .le. 85  .and. &
                loni(k) .le. 330 .and. &
                loni(k) .gt. 285      ) then
               do l = 1, nsoilay
                  tsoi(k,l) = xint(k,1) + 273.16
                  wisoi(k,l) = 1.0
               end do
               tg(k) = tsoi(k,1)
               hsno(k,1) = 0.05
               hsno(k,2) = 2.
               hsno(k,3) = 2.
               fi(k) = 1.
            end if
    
            if (lati(k) .ge. 66  .and.  &
                lati(k) .le. 85  .and.  &
                loni(k) .le. 345 .and.  &
                loni(k) .gt. 330) then
               do l = 1, nsoilay
                  tsoi(k,l) = xint(k,1) + 273.16
                  wisoi(k,l) = 1.0
               end do
               tg(k) = tsoi(k,1)
               hsno(k,1) = 0.05
               hsno(k,2) = 2.
               hsno(k,3) = 2.
               fi(k) = 1.
            end if
         end do
      end if
      return
end subroutine coldstart
