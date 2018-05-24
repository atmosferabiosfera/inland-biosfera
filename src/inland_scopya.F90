#include "inland_config.h"
!----------------------------------------------------------------------
! copies array arr to brr,for 1st nt words of arr
subroutine scopya (nt, arr, brr)
!----------------------------------------------------------------------
      implicit none

! Arguments
      integer nt     
      real*8 arr(nt), &  ! input
             brr(nt)     ! output

! Local variables
      integer ia

      do 100 ia = 1, nt
         brr(ia) = arr(ia)
100   continue
      return
end subroutine scopya
