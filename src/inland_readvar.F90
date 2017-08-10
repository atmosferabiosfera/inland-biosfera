#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! .....................................................................
! readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
!.....................................................................

subroutine readvar(filen,varname,values,istart,icount,nlpt,ierror)
!
      use inland_comwork, only: cdummy,nlonsub,nlatsub
      use inland_control, only: env_debug
!
      implicit none
!
! INPUT
! filen - character*(*) - file name from which to read data
! varname - character*(*) - name of variable from which to read
! Ignored if varname variable is only 3-d
! istart - integer(*) - starting points along each dimension of
! the variable, for example the 1st point a 4-d variable is (1,1,1,1)
! and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
! icount - integer(*) - number of points along each dimension to read,
! for example, to read in a single lat/lon grid from a 4-d variable,
! icount would be (nlon,nlat,1,1)
! nlpt - number of tiles to copy 
!  -1 do not copy anything (must call arr2vec afterwards)
!   0 replicate to all tiles using arr2vec
!  >0 replicate to nlpt tiles using arr2vec_tile 
!
! OUTPUT
! values - real(*) - returned real values of the designated hyperslab
! ierror - integer - error code = 0 if no error, < 0 if error occurred
! note that dimension values (alons,alats,vals3d,times)
! must now be read seperately, see readit()
      character*(*) filen, varname
      integer :: istart(*), icount(*), ierror, nlpt
      real*8 :: values(*)

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

      integer :: idies, idvar, ndims, idlon, idlat, id3d, ierr, itype
      integer :: idtime, ilpt

      ierror = 0

! Open file, but don't fail if already open
! -----------------------------------------
      ierr = NF_OPEN(filen,NF_NOWRITE,idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in readvar '//trim(filen)
         print *, NF_STRERROR(ierr)
      end if

! Get id of variable
! ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 1 in readvar'
         print *, NF_STRERROR(ierr)
         goto 999
      end if

! Inquire number of dimensions
! ----------------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 2 in readvar'
         print *, NF_STRERROR(ierr)
         goto 999
      end if

! Get variable
! ------------
      if ( env_debug .gt. 1 ) print *,'readvar '//trim(varname)//' from '//trim(filen)
      if ( nlpt .gt. -1 ) then
         ierr = NF_GET_VARA_DOUBLE(idies,idvar,istart,icount,cdummy)
      else
         ierr = NF_GET_VARA_DOUBLE(idies,idvar,istart,icount,values)
      end if
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 3 in readvar '//varname
         print *, NF_STRERROR(ierr)//' (#',ierr,')'
         ! special case - if error is "Index exceeds dimension bound" (-40) print error but do not exit
         ! this is a work-around for reading tilefrac with no waterprop in readit
         if ( ierr .eq. NF_EINVALCOORDS ) then
            print *, 'this error may not be fatal...'
            ierror = -1
            goto 9999
         else
            goto 999
         end if
      end if

! Get values of dimension variables
! note that dimension values (alons,alats,vals3d,times)
! must now be read seperately, see readit()
! ---------------------------------

! Close file
! ----------
      ierr = NF_CLOSE(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 11 in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if


! Copy values to each tile if requested
! -------------------------------------
      if ( nlpt .gt. 0 ) then
         do ilpt = 1, nlpt
            call arr2vec_tile(cdummy((ilpt-1)*nlonsub*nlatsub + 1), values, ilpt)
         end do
      else if ( nlpt .eq. 0 ) then
         call arr2vec(cdummy, values)
      end if ! nlpt

!     return without error
!     --------------------
     return

!     Got error - exit
!     ----------------
999   write(*,*) 'varname = '//varname
      stop 1
9999  write(*,*) 'varname = '//varname

end subroutine readvar
