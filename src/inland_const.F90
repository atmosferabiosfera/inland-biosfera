#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine const (arr, nar, value)
! ---------------------------------------------------------------------
      implicit none
!
! sets all elements of real vector arr to value
!
      integer j, nar
      real*8 arr(nar), value
!
      do 100 j = 1, nar
        arr(j) = value
100   continue
!
      return
end subroutine const
