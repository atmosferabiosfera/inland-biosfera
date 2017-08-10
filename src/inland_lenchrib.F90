#include "inland_config.h"
integer function lenchrib(chrstg)
!-----------------------------------------------------------------------
!
! Return position of right-most non-blank, non-null character
! in chrstg.
!
!---------------------------Code history--------------------------------
!
! Original version:  L. Bath, April 1992
! Standardized:      L. Bath, June 1992
!                    J. Rosinski April 1994
!                    T. Acker, March 1996
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
      character*(*) chrstg       !  Input character string
!
!--------------------------Local Variables------------------------------
!
      integer l     ! loop counter
!
!-----------------------------------------------------------------------
!
      lenchrib = 0
      do l=len(chrstg),1,-1
!
! Find right-most non-blank character in string
!
         if (chrstg(l:l).ne.' ' .and. chrstg(l:l).ne.char(0)) then
            lenchrib = l
            goto 10
         end if
      end do
!
10    return
!
end function lenchrib


