#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ies-io.f last update 9/9/98 Christine C. Molling, cmolling@facstaff.wisc.edu
! University of Wisconsin - Madison; Institute for Environmental Studies;
! Center for Climatic Research; Climate, People, and Environment Program
!
! You are free to change or distribute my code as long as you
! 1) acknowlege me, Christine C. Molling (cmolling@facstaff.wisc.edu)
! as the author of the original code and 
! 2) do not sell it.
! Of course if there is anything wrong with the code, or you somehow
! encounter damage to your reputation/wallet because of its use, you
! are forbidden to sue me or anyone else.
!
!     These subroutines can be used as a primary interface to the ies
! format netcdf files.  You can also use the slightly lower level functions
! found in ies.f or the low level netcdf commands (see the netcdf V3
! users guide).
!SUBROUTINE LIST:
! inifile - create a new file with dimensions and variables for
!  longitude, latitude, an optional 3rd real dimension, and time
!  Time must have units of days since a date!!!
! inifilec - create a new file with dimensions and variables for
!  longitude, latitude, a 3rd character-variable dimension, and time
!  Time must have units of days since a date!!!
! inivar - initialize a new variable
! endini - end initialization of file: sync file to disk, close file
! writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
! readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
!.....................................................................
! writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
!.....................................................................

! subroutine which takes care of basic netcdf i/o, without tiles
subroutine writevar_notile(filen,idies,varname,values,istart,icount,times,ierror)
!
      implicit none
!
#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif
!
!INPUT
!     filen - character*(*) - file name to which to write data
!     varname - character*(*) - name of variable to which to write
!     istart - integer(*) - starting points along each dimension of
!      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
!      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
!     icount - integer(*) - number of points along each dimension to write,
!      for example, to write a single lat/lon grid to a 4-d variable,
!      icount would be (nlon,nlat,1,1).
!     values - real(*) - real values of the designated hyperslab
!     times - real(*) - time vector vector for those points written (ignored
!      if variable is 2-d).
!OUTPUT
!     ierror - integer - error code = 0 if no error, < 0 if error occurred

!FIXME: These asterisks on var dimensions should be : to comform to fortran 90!
      character*(*) filen, varname
      integer istart(*), icount(*), ierror
      real*8 values(*), times(*)

      integer idies, idvar, ndims, ierr, idtime, idtw, iddate

      ierror = 0

!     Get id of variable
!     ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 1 in writevar_notile ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if

if ( trim(varname) .eq. "zbot") print *,"idvar:",idvar

!     Inquire number of dimensions
!     ----------------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 2 in writevar_notile ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if
         
!     Put variable
!     ------------
      ierr = NF_PUT_VARA_DOUBLE(idies,idvar,istart,icount,values)

      if (ierr .ne. NF_NOERR) then
         print *, 'Error 3 in writevar_notile ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if

!     Put time value(s) (removed time weights and dates - useless)
!     --------------------------------------------------
      if (ndims .gt. 2) then
         ierr = NF_INQ_VARID(idies,'time',idtime)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error 4 in writevar_notile ',varname
            print *, NF_STRERROR(ierr)
            ierror = -1
            goto 999
         end if
         ierr = NF_PUT_VARA_DOUBLE(idies,idtime,istart(ndims),icount(ndims), &
                                   times)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error 5 in writevar_notile ',varname
            print *, NF_STRERROR(ierr)
            ierror = -1
            goto 999
         end if
      end if

!     return without error
!     --------------------
      return

!     Got error - exit
!     ----------------
999   write(*,*) 'varname = '//varname
      stop 1

end subroutine writevar_notile


! subroutine which takes care of tiles, including last tile (average)
subroutine writevar(filen,idies,varname,values,istart,icount,times,ierror)

      use inland_comwork, only: cdummy
      use inland_parameters
      use inland_subgrid

      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

!INPUT
      character*(*) filen, varname
      integer istart(*), icount(*)
      real*8 values(*), times(*)
!OUTPUT
      integer ierror

      integer idies, idvar, ndims, ierr, idtime, idtw, iddate
      integer i, k, k2, ilpt, nlpt, ixdim, nxdim

      real*8  nodata
      real*4  nodata2
      integer xtype, nodataint
      character*80 chtmp
      integer dimids(5)

!     Get id of variable
!     ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 1 in writevar ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if

if ( trim(varname) .eq. "zbot") print *,"idvar:",idvar

!     Get dimensions info
!     -------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 2 in writevar ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if
      ierr = NF_INQ_VARDIMID(idies,idvar,dimids)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 3 in writevar ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if

      ! get tile and extra dimension count (pft/nsoilay/etc)
      ierr = NF_INQ_DIMNAME(idies,dimids(ndims-1),chtmp)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 4 in writevar ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         goto 999
      end if
      nxdim = 1
      ! check if we have a tile dimension
      if ( chtmp .ne. 'tile' ) then
         nlpt = 1
         if ( ndims .eq. 4 ) nxdim = icount(ndims-1)
      else
         nlpt = mlpt
         if ( ndims .eq. 5 ) nxdim = icount(ndims-2)
      end if

if ( trim(varname) .eq. "zbot") print *,"dims:",nxdim

!     prepare cdummy array
!     -----------------------------------------           
      nodata=ocean

      ! get missing_value for this variable, assume ocean on error
      ierr = NF_INQ_VARTYPE(idies,idvar,xtype)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning 2 in writevar '//trim(varname)
         print *, NF_STRERROR(ierr)
      else
         if ( (xtype .eq. NF_INT) .or. (xtype .eq. NF_BYTE) ) then
            ierr = NF_GET_ATT_INT(idies,idvar,'missing_value',nodataint)
            nodata=nodataint
         else if ( (xtype .eq. NF_FLOAT) ) then
            ! ugly hack to get nodata to work when type is float
            ierr = NF_GET_ATT_REAL(idies,idvar,'missing_value',nodata2)
            nodata = nodata2
         else
            ierr = NF_GET_ATT_DOUBLE(idies,idvar,'missing_value',nodata)
         end if
         if (ierr .ne. NF_NOERR) then
            nodata=ocean
            print *, 'Warning 3 in writevar '//trim(varname)
            print *, NF_STRERROR(ierr)
         end if
      end if

      ! don't fill with nodata anymore - this is too costly when writing to daily/monthly files

      ! do a simple copy if there are no tiles and extra dims
      if ( ( ndims .lt. 4 ) .or. &
           ( ( ndims .eq. 4 ) .and. ( icount(ndims-1) .eq. 1 ) ) ) then
         
         call vec2arr(values, cdummy(1), 1, nodata)
    
      ! else copy all tiles 
      else 

         ! copy each tile to cdummy
         ! TODO this is hard-coded to mlpt, instead of icount(ndims-1)
         ! this not great, but simpler because of average
         do ilpt = 1, mlpt
            do ixdim = 1, nxdim
	if ( trim(varname) .eq. "zbot") print *,ilpt,ixdim
               call vec2arr2d(values, &
                              cdummy((ilpt-1)*nxdim*nlonsub*nlatsub + &
                              (ixdim-1)*nlonsub*nlatsub+1), &
                              ilpt, nodata,ixdim,nxdim)
            end do ! ilpt
         end do ! ixdim

         ! calculate average if third dim has mlpt+1 length
         if ( ( mlpt .gt. 1 ) .and. ( icount(ndims-1) .eq. mlpt+1 ) ) then
            
            do ixdim = 1, nxdim

               buffer1(:)=0.
      
               do ilpt = 1, mlpt
                  ! need to loop over every point to calculate average
                  do i=1,npoi1
                     k=subgrid_get_index(i,ilpt)
                     if ( k .eq. 0 ) then
                        write (*,*) 'ERROR in writevar(), got subgrid index',k,'with i=',i,' ilpt=',ilpt
                     else
                        k2 = (ixdim-1)*mlpt*npoi1 + k
                        ! TODO check for nodata, must test more
                        if ( ( values(k2) .ne. nodata ) .and. ( tilefrac(k) .gt. 0.0 ) ) then
                           ! special case for tilefrac and burnarea - write total instead of average
                           ! some points the total is 0.9999* but no big deal
                           if ( varname .eq. "tilefrac" .or. varname .eq. "burnarea" ) then
                              buffer1(i) = buffer1(i) + values(k2)
                           else
                              buffer1(i) = buffer1(i) + values(k2) * tilefrac(k)
                           end if
                        end if
                     end if
                  end do ! i
               end do ! mlpt
               
               ! copy average to last tile
               call vec2arr(buffer1, &
                            cdummy((mlpt)*nxdim*nlonsub*nlatsub + &
                            (ixdim-1)*nlonsub*nlatsub + 1), &
                            1, nodata)
               
            end do !ixdim          
            
         end if ! calculate average

      end if ! no tiles/dims

!     call writevar with cdummy
!     -------------------------
      call writevar_notile(filen,idies,varname,cdummy,istart,icount,times,ierror)
      if (ierror .ne. 0) goto 999

!     return without error
!     --------------------
      return

!     Got error - exit
!     ----------------
999   write(*,*) 'varname = '//varname
      stop 1

end subroutine writevar
