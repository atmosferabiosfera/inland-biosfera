#include "inland_config.h"
! ----------------------------------------------------------------------------
subroutine rotation(irestart, iyrrestart)
! ----------------------------------------------------------------------------
!
! subroutine used to plant different crop types in different years
! to simulate typical crop rotations in the United States 
! could be modified in future if rotations would include natural
! vegetation types such as grasslands or biofuel (switchgrass) crops
!
! currently three main types of rotations are used
! if iwheattype eq:
! 2: maize/soybean rotation (alternating)
! 3: maize/soybean/spring wheat
! 4: soybean/winter wheat/maize
! note: by doing continuous winter wheat, land is fallow 
! only from harvest (june-july) through planting (Sept-Nov).
!
!
! common blocks
!
      use inland_parameters
      use inland_control,only:iyear
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr

!
! local variables
!
      real*8 xfrac
!
      integer irestart, &
              iyrrestart, &
              iyrrot,   &
              i,        &
              idiv
!
! begin grid
! in this case, irotation also is the number of crops in a specified rotation
!
! assumes that natural vegetation existence arrays are set to 0.0
! in climate.f (existence) 
!
      if (irestart .eq. 1) then
         iyrrot = iyrrestart
      else 
         iyrrot = 1950  ! base year to start crop rotations
      endif
!
! look at the fraction remainer to determine which crop should
! be planted
!
      if (irotation .eq. 1) idiv = 3
      if (irotation .eq. 2) idiv = 2
      if (irotation .eq. 3) idiv = 3
      if (irotation .eq. 4) idiv = 3
!
      xfrac = mod ((iyear - iyrrot), idiv)  

      do 100 i = lbeg, lend
!
! 2:
! two-crop rotation (standard soybean/corn)
! alternate between even / odd years
!
         if (irotation .eq. 2) then
            if (xfrac .eq. 0) then
               isoybean = 1
               imaize   = 0
               iwheattype   = 0
               exist(i,13) = 1.0
               exist(i,14) = 0.0
               exist(i,15) = 0.0
            else
               isoybean = 0
               imaize  = 1
               iwheattype  = 0
               exist(i,13) = 0.0
               exist(i,14) = 1.0
               exist(i,15) = 0.0
            endif
! 3:
! rotation with 3 crops (corn, soybean, spring wheat) 
!
         elseif (irotation .eq. 3) then
            if (xfrac .eq. 0) then 
               isoybean = 0
               imaize  = 1
               iwheattype   = 0
               exist(i,13) = 0.0
               exist(i,14) = 1.0
               exist(i,15) = 0.0
            elseif (xfrac .eq. 1.0) then
               isoybean = 1
               imaize = 0
               iwheattype  = 0
               exist(i,13) = 1.0
               exist(i,14) = 0.0
               exist(i,15) = 0.0
            else
               isoybean = 0
               imaize  = 0
               iwheattype  = 1
               exist(i,13) = 0.0
               exist(i,14) = 0.0
               exist(i,15) = 1.0
            endif

! 4:
! 3 crop rotation with winter wheat and soybean planted in same year
! winter wheat harvested in year 2
! maize grown in year 3
 
         elseif (irotation .eq. 4) then 

! soybean planted/harvested  
! winter wheat planted

            if (xfrac .eq. 0.0) then 
               isoybean = 1
               imaize  = 0
               iwheattype  = 2
               exist(i,13) = 1.0
               exist(i,14) = 0.0
               exist(i,15) = 1.0

! winter wheat harvested in year 2

            elseif (xfrac .eq. 1.0) then
               isoybean = 0
               imaize = 0
               iwheattype   = 2
               exist(i,13) = 0.0
               exist(i,14) = 0.0
               exist(i,15) = 1.0
     
! maize planted/harvested in year 3

            else
               isoybean = 0
               imaize   = 1
               iwheattype  = 0
               exist(i,13) = 0.0
               exist(i,14) = 1.0
               exist(i,15) = 0.0
            endif
         endif
100   continue
!
      return
!      
end subroutine rotation

