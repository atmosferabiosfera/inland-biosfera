#include "inland_config.h"
module inland_comhour

      implicit none
      public
      save

!
! --------
! comhour.h  include file  last update 10.02.02 C. Kucharik 
! --------
!
! this include file is to hold values for variables related to inputting
! a data file of meteorological quantities (hourly) to be used to drive
! IBIS using real weather data
!
      integer ndpts, imetyear, dmetyear, imetend, dmetend
!
      parameter (ndpts=9000) ! the number of data points in the file 
!
!
      real*8, dimension(:), allocatable :: var1, var2, var3, var4, var5, var6
!             var1(ndpts),
!             var2(ndpts),
!             var3(ndpts),
!             var4(ndpts),
!             var5(ndpts),
!             var6(ndpts)

!
end module inland_comhour
