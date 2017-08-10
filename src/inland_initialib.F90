#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine initialib (isimveg, irestart, iyearp, imonthp, rdlsf, lati,loni)
! ---------------------------------------------------------------------
!     intializes the INLAND variables or reads the restart files (at
!     the beginning of the year or month)
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      logical rdlsf        ! true: initial data for Antarctica, Greenland.
      integer isimveg      ! 0: static veg, 
                           ! 1: dynamic veg, 2: dynamic veg-cold start
      integer irestart     ! 0: not restart of INLAND, 1 restart of INLAND
      integer iyearp       ! initial year of simulation 
                           ! (equals iyrlast if restart)
      integer imonthp      ! initial month of simulation 
                           ! (equals imolast if restart)

      real*8 lati(npoi)       ! latitude of land point on long vector [rads]
      real*8 loni(npoi)       ! longitude of land point on long vector [rads]

      if (irestart .eq. 0) then
         call coldstart (rdlsf,  lati, loni)
      else
#ifdef SINGLE_POINT_MODEL
! 0D model version does NOT support restart!
         print *,"0D version does NOT support restart."
         stop 1
#else
         !TODO: test restart after the model runs correctly from the beginning
         ! Just work with isimveg 1, with dynamic vegetation
!this part is for calculating the crop 0D
    if(npoi .gt. 1)then
         if(isimveg.ge.1) call restartib (iyearp, imonthp)
    endif
#endif /* SINGLE_POINT_MODEL */
      end if

! initialize physical consts, dimensions, unit numbers, lsx model
      call inisurf (irestart)

! initialize snow model
      call inisnow 

! initialize soil model
      call inisoil (irestart)

! initialize vegetation parameters
      call iniveg (isimveg, irestart)

! initialize variables for time averaging
      call inisum (irestart)
      return
end subroutine initialib
