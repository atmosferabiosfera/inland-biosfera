#include "inland_config.h"
! ---------------------------------------------------------------------
      subroutine constint (arr, nar, value)
! ---------------------------------------------------------------------
!
! sets all elements of real vector arr to value
!
!
! Arguments
!
      integer nar, value
      integer arr(nar)
!
! Local variables
!
      integer j
!
      do 100 j = 1, nar
        arr(j) = value
 100  continue
!
      return
end subroutine constint
