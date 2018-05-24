#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine commenthandler (parm_unit, parm_file)     
!-----------------------------------------------------------------------
! #, *, c, C, !, // denote comments in the parameter file. Therefore the 
! input reading routine must ignore these symbols and any text to their right. 
! This subroutine checks for the occurrence of these symbols, whenever RD_PARAM
! detects an error reading numeric data. If the error was caused by text starting
! with one of these symbols, the text is skipped over and the routine returns
! so that the next line of data can be read in.

! If an error is detected, the line containing the offending characters is 
! echoed to the screen and execution stops. But if someone were so-inclined,
! they could write some extra code to analyse the problem and possibly recover
! from it. 
!---------------------------Code history--------------------------------
!
!     Original version: David Price
!
! ---------------------------------------------------------------------
      implicit none

! Parameters
      integer*4 parm_unit
      character*(*)  parm_file ! name of file from which input is being read in

! Local variables
      integer*4 iocode   ! dummy variable required by system for error handling?
      integer*2 is_comment

! arbitrary string. Input lines cannot exceed 255 chars!
      character*255 inputline 

      backspace(parm_unit) ! back up to start of field that caused this error
 
      read (parm_unit, 1002, end=800, err=801, iostat=iocode) inputline
1002  format (A255)

      is_comment = 0          ! must initialize this correctly

      if ((index (inputline, 'C') .eq. 1) .or. &
          (index (inputline, 'c') .eq. 1) .or. &
          (index (inputline, '*') .eq. 1) .or. &
          (index (inputline, '#') .eq. 1) .or. &
          (index (inputline, '!') .eq. 1) .or. &
          (index (inputline, '/') .eq. 1)) is_comment = 1

      if (is_comment .NE. 1) goto 801 ! we have a problem....
 
      return  ! But if we get to here, all's well and we can move on.

! FIXME: the following three commands are never reached in life
800   write (*, 1000) parm_file, inputline ! For now, just echo the line and stop
1000  format ('CommentHandler: Unexpected EOF encountered in ', &
              A12, '\nLast input was:', A255)
      stop 1

801   write (*, 1001) parm_file, inputline ! For now, just echo the line and stop.
1001  format ('CommentHandler: Error encountered reading data ', &
              'from ', A12, '\nLast input was:', A255)
      stop 2
end subroutine commenthandler
