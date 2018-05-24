#include "inland_config.h"
! ---------------------------------------------------------------------------------
! added for Green-Ampt infiltration 
!
       subroutine deltasw(accuminf, dtsw, layerno, i)
!
! ---- calculate change in volumetric water content from surface to wetting front          
!
      use inland_parameters
       use inland_comsoi
       use inland_comatm
!
       real*8 accuminf,	   & !actual accumlative infiltration (mm)
     	    dtsw,  	   & ! weighted change in volumetric soil water content corresponding to wetting depth 
     	    temp,  	   & ! temperary variable
     	    depth, 	   & ! wetting depth in the soil layer in which wetting front is located
     	    wetdep, 	   & ! wetting depth from surface to wetting front
	    dt_sw(nsoilay)   !delta volumetric content for six soil layer
!
       integer layerno,  & ! Layer number
       	       num,   	 & ! total number of layers from top layer to the layer of wetting front
     	       i      	   ! grid cell number
!       
       do j = 1, nsoilay 
              dt_sw(j) = poros(i, j) * (1.0 - sice(i, j)) * (1.0 - swater(i, j))  !original-ga
		if (j .eq. 1) then
		endif
             dt_sw(j) = max(0.0000001, dt_sw(j)) 
       end do
!
       temp = 0.0     ! initialization
       wetdep = 0.0   ! initialization
! 
       do j = 1, nsoilay 
          temp = temp + dt_sw(j) * hsoi(j) * 1000.0
          num = j             
          if (temp .ge. accuminf) goto 100 
       end do
!
100    depth = (accuminf - (temp - dt_sw(num) * hsoi(num) * 1000.0))/(dt_sw(num) * 1000.0)
       depth = max(0.0, depth)
       temp = 0  ! renitialization    
!       
       do j = 1, num - 1
          wetdep = wetdep + hsoi(j)     
       end do
!
       wetdep = wetdep + depth
!
       do j = 1, num - 1
          temp = temp + dt_sw(j) * (hsoi(j) / max(0.0000001, wetdep))
       end do 
       dtsw  = temp + dt_sw(num) * (depth / max(0.0000001, wetdep))
       layerno = max(num, 1)
       layerno = min(num, nsoilay)
       dtsw = max(0.0000001, dtsw)
       return
       end subroutine deltasw
