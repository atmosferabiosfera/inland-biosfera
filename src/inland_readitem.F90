#include "inland_config.h"
!*******************************************************************************

subroutine readitem (funit, fname, item)
!-----------------------------------------------------------------------
! Simple routine to read in data in free format, with a built-in error handler 
! designed to locate and skip over comment strings buried in the data stream.
! ---------------------------------------------------------------------
      implicit none
      integer*4 funit,  &
                iocode   ! dummy variable required by system for error handling?

      real*8    item
      character*(*) fname

101   read (funit, *, end=999, err=911, iostat=iocode) item
      return

911   call commenthandler (funit, fname)
      goto 101

999   write (*,2) fname, funit
2     format ('RD_PARAM: Unexpected EOF encountered in file ',A20,' on unit ', &
              I2)
      stop
end subroutine ReadItem
