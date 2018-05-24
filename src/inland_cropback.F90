#include "inland_config.h"
! ---------------------------------------------------------------------
      subroutine crop_back (arr3,arr2,j,vrn)
! ---------------------------------------------------------------------
!
! chooses between two things.  Used in inland_canopy.F90
!
!
      use inland_parameters

      real*8 arr3(lbeg:lend,npft,60), arr2(lbeg:lend,npft)

      integer  i,vrn,j
!
      do i = lbeg, lend
         arr2(i,j) = arr3(i,j,vrn)
      enddo
!
      return
      end subroutine crop_back
