#include "inland_config.h"
!----------------------------------------------------------------------
! FIXME: useless subroutine to fill an array with a given value.
subroutine setarrib (pa, kdim, pvalue)
!----------------------------------------------------------------------

      implicit none

! ------------------------ code history ---------------------------
! source file:       setarrib.F
! purpose:           set array pa(kdim) to pvalue
! date last revised: March 1996 - lsm version 1
! author:            Gordon Bonan
! standardized:      J. Truesdale, Feb. 1996
! reviewed:          G. Bonan, Feb. 1996
! -----------------------------------------------------------------

! ------------------------ input/output variables -----------------
! input
      integer kdim      !dimension of array pa
      real*8 pvalue       !value to store in pa

! output 
      real*8 pa(kdim)     !array to set
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
      integer j         !loop index
! -----------------------------------------------------------------
      do j = 1, kdim
         pa(j) = pvalue
      end do
! 
      return
!      
end subroutine setarrib
