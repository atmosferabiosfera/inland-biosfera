#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
      subroutine arr3_arr2 (arr3,arr2,lre)
! ---------------------------------------------------------------------
!
! chooses between two things.  Used in inland_canopy.F90
!
!
      use inland_parameters
	
      real*8 arr3(lbeg:lend,2100), arr2(lbeg:lend)
	integer  i,lre
!
	do i=lbeg, lend
	 arr2(i) = arr3(i,lre)
	enddo
!
      return
      end subroutine arr3_arr2
