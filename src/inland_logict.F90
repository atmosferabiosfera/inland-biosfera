#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine logict (arrl, nar)
! ---------------------------------------------------------------------
      implicit none
!
      integer j, nar
      logical arrl(nar)
!
      do 100 j = 1, nar
         arrl(j) = .TRUE.
 100  continue
      return
end subroutine logict
