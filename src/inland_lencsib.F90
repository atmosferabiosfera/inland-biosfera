#include "inland_config.h"
integer function lencsib (chrstr)
! ------------------------ code history ---------------------------
! source file:       lencsib.F
! purpose:           return position of right-most non-blank, non-null
!                    character in chrstr
! date last revised: March 1996 - lsm version 1
! author:            Gordon Bonan
! standardized:      J. Truesdale, Feb. 1996
! reviewed:          G. Bonan, Feb. 1996
! -----------------------------------------------------------------
      implicit none
! ------------------------ input/output variables -----------------
! input
      character*(*) chrstr       !input character string
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
      integer :: l
! -----------------------------------------------------------------

      lencsib = 0
      do 100 l = len(chrstr),1,-1
         if (chrstr(l:l).ne.' ' .and. chrstr(l:l).ne.char(0)) then
            lencsib = l
            goto 10
         end if
100   continue
10    return
end function lencsib
