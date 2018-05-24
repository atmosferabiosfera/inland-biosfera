#include "inland_config.h"
! Kaiyuan Li for Green-Ampt infiltration NOTE NEW FUNCTION
! caluclate change in volumetric water content from surface to wetting front          

! ---------------------------------------------------------------------
subroutine delta_sw(accuminf, dtsw, layerno, i)
! ---------------------------------------------------------------------

      use inland_parameters
      use inland_comatm     ! Castanho Kai NOTE var not used
      use inland_comsoi
      use inland_comveg     ! Castanho Kai NOTE var not used   

      implicit none

! input variables
! accuminf,   !actual accumlative infiltration (mm)
! dtsw,       ! weighted change in volumetric soil water content corresponding to wetting depth 

       real*8 accuminf, dtsw

!      layerno,    ! Layer number
!      i           ! grid cell number
       
       integer layerno, i, j

! local variables

       real*8 dt_sw(nsoilay)  ! delta volumetirc content for six soil layer
       real*8 temp, &         ! temperary variable
              depth, &        ! wetting depth in the soil layer in which wetting front is located
	      wetdep          ! wetting depth from surface to wetting front
	      
       integer  num           ! total number of layers from top layer to the layer of wetting front

       do j = 1, nsoilay 
          dt_sw(j) = poros(i, j) * max(dble(0.0), (1.0 - sice(i, j))) * max(dble(0.0), (1.0 - swater(i, j)))
       end do

       temp = 0.0     ! initializaton
       wetdep = 0.0   ! initialization
 
       do j = 1, nsoilay 
          temp = temp + dt_sw(j) * hsoi(j) * 1000.0
          num = j
          if (temp .ge. accuminf) goto 100 
       end do

100    depth = (accuminf - (temp - dt_sw(num) * hsoi(num) * 1000.0))/(dt_sw(num) * 1000.0)
       depth = max(dble(0.0), depth)
   
       temp = 0  ! renitialization    
      
       do j = 1, num - 1
          wetdep = wetdep + hsoi(j)     
       end do

       wetdep = wetdep + depth
       do j = 1, num - 1
          temp = temp + dt_sw(j) * (hsoi(j) / max(0.0000001, wetdep))
       end do 
       dtsw  = temp + dt_sw(num) * (depth / max(0.0000001, wetdep))
       layerno = max(num, 1)
       layerno = min(num, nsoilay)
       dtsw = max(dble(0.0), dtsw)
   
end subroutine delta_sw    
