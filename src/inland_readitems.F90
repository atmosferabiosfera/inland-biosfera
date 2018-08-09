#include "inland_config.h"
!*******************************************************************************
subroutine readitems (funit, fname, N, item1, item2, item3, item4, &
                      item5, item6, item7, item8, item9, item10, item11, item12)
!-----------------------------------------------------------------------
! Hokey routine to read in up to N data items in free format, with a built-in 
! error-handler designed to locate and skip over comment strings buried in the
! data stream. Item1 to Item10 are holders for a sequence of legitimate data
! values.  N specifies how many of these are to be used. Therefore N must be 
! less than or equal to 10 (unless you increase the number of items). Note
! that N = 1 is allowed, but it is easier to use subroutine ReadItem for 
! reading single values. 
! ---------------------------------------------------------------------

      implicit none
     
      integer*4 funit,  &
                iocode, &   ! dummy variable required by system for error handling?
                N

      real*8 item1, item2, item3, item4, item5, item6, item7, item8, item9, item10, item11, item12
      character*(*) fname

! Check to see that N is within acceptable range. Stop if not.
      if ((N.gt.12).or.(N.lt.0)) then
         write (*,3)
3        format ('READITEMS: Invalid number of items specified.')
         stop
      end if

100   goto (101,102,103,104,105,106,107,108,109,110,111,112) N
  
101   read (funit, *, end=999, err=911, iostat=iocode) item1
      return

102   read (funit, *, end=999, err=911, iostat=iocode) item1, item2
      return

103   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3
      return

104   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4
      return

105   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5
      return

106   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6
      return

107   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6, item7
      return

108   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6, item7, item8
      return

109   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6, item7, item8, item9
      return

110   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6, item7, item8, item9, item10
      return
      
111   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6, item7, item8, item9, item10, item11
      return

112   read (funit, *, end=999, err=911, iostat=iocode) item1, item2, item3, item4, item5, item6, item7, item8, item9, item10, item11, item12
      return

911   call commenthandler (funit, fname)
      goto 100

999   write  (*,2) fname, funit
2     format ('RD_PARAM: Unexpected EOF encountered in file ', A20, ' on unit ', I2)
      stop
end subroutine readitems
