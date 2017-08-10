#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine logicf (arrl, nar)
! ---------------------------------------------------------------------
      implicit none
!
      integer j, nar
      logical arrl(nar)
!
      do 100 j = 1, nar
         arrl(j) = .FALSE.
100   continue
!
      return
end subroutine logicf
