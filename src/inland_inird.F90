#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! This subroutine obtains the year+1 of the year in the units attribute
! of a file.  Low-level netcdf commands are used.
! ---------------------------------------------------------------------
subroutine inird(file,istyr)
! ---------------------------------------------------------------------
      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

! Arguments
      character*(*) file
      integer istyr       

! Local Variables
      integer idies, istat, idtime, lf1

! Remove this block for Absoft f77, D.Pol 15-Mar-2002
! unsure if needed for SGI version
!      integer NF_OPEN,        ! netcdf function
!     >        NF_NOWRITE,     ! '
!     >        NF_NOERR,       ! '
!     >        NF_INQ_VARID,   ! '
!     >        NF_GET_ATT_TEXT ! '
!              
      character*80 units
! ---------------------------------------------------------------------
! open file
      istat = NF_OPEN(file,NF_NOWRITE,idies)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to open file'
         print *, file
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if

! get units attribute for time
      istat = NF_INQ_VARID(idies,'time',idtime)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to get time id'
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
      units = ' '
      istat = NF_GET_ATT_TEXT(idies,idtime,'units',units)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to get time units'
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if

! put character units year into integer variable, add 1
      lf1 = index(units,'since') + 6
      read(units(lf1:lf1+3),'(i4)') istyr
!      read(units(12:15),'(i4)') istyr
!      istyr = istyr + 1
      return
end subroutine inird
