#include "inland_config.h"
subroutine co2 (co2init)
! ---------------------------------------------------------------------
      use inland_comatm, only: co2conc
      use inland_control, only: iyear
      implicit none

! Arguments 
      real*8 co2init !initial co2 concentration in mol/mol (coming from inland.infile)

! Variables
      integer i      !number of dimensions of variable 'dummy'
      real*8, dimension(2) :: dummy !temporary variable to hold data
                                    !read from input file 'co2.data.txt'

! First, assign co2 concentration as defined in inland.infile (co2init)
      co2conc = co2init

! Then, read co2 concentration for this year. Nothe that all year and co2
! concentration values in the file will be read, but the co2conc will be 
! assigned only when current simulation year is equal to year read in the file
      open (unit=20,file='input/co2.data.txt',status='old')
15    read (20,*,end=101)(dummy(i),i=1,2)
      if (dummy(1).eq.iyear) then
         co2conc = dummy(2)
      endif
      goto 15
101   close (20)

! 1992 IPCC estimates
!       iyr = iyear - 1860 + 1
!       co2conc = (297.12 - 0.26716 * iyr +
!    >                      0.0015368 * iyr**2 +
!    >                      3.451e-5 * iyr**3) * 1.e-6
!
! M. El Maayar: 1996 IPCC estimates
!       iyr = iyear - 1860 + 1
!       co2conc = (303.514 - 0.57881 * iyr +
!    >                      0.00622 * iyr**2 +
!    >                      1.3e-5 * iyr**3) * 1.e-6
      return
end subroutine co2
