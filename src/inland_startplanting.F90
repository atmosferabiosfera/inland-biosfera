#include "inland_config.h"
! ----------------------------------------------------------------------------
subroutine startplanting (iyear0,iyear)
! ----------------------------------------------------------------------------
!
! common blocks
!
      use inland_ibisparm
      use inland_comcrop
!
      integer i,iyear0,iyear       
!
      if (nratoon.eq.4) then
         do i=1,npoi
            if(mod(i,6) .eq. 0 .and. iyear .ge. iyear0)     ncyears(i) = ncyears(i) +1 
            if(mod(i,6) .eq. 1 .and. iyear .ge. iyear0 + 1) ncyears(i) = ncyears(i) +1 
            if(mod(i,6) .eq. 2 .and. iyear .ge. iyear0 + 2) ncyears(i) = ncyears(i) +1 
            if(mod(i,6) .eq. 3 .and. iyear .ge. iyear0 + 3) ncyears(i) = ncyears(i) +1 
            if(mod(i,6) .eq. 4 .and. iyear .ge. iyear0 + 4) ncyears(i) = ncyears(i) +1 
            if(mod(i,6) .eq. 5 .and. iyear .ge. iyear0 + 5) ncyears(i) = ncyears(i) +1 
         enddo
     endif
!
     return
! 
end subroutine startplanting
