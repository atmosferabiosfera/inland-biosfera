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
!
!.....................................................................
! inifile - create a new file with dimensions and variables for
!  longitude, latitude, a number of extra dimensions, and time.
!  Time must have units of days since a date!!!
!.....................................................................

! version which initializes a file with numXdim extra dimensions
subroutine inifile_dims(idies,filen,title,source,history,nlon,alons,nlat,alats, &
                    nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXvals, &
                    tunits,calendar,ierror)
!
       use inland_control, only: env_compressout
!
!INPUT
!     filen - character*(*) - name for new file
!     title - character*(*) - title for file
!     source - character*(*) - source for file
!     history - character*(*) - date of creation, and possibly author
!     nlon - integer - number of point along longitude direction
!     alons - real(nlon) - values of longitudes
!     nlat - integer - number of point along latitude direction
!     alats - real(nlat) - values of latitudes
!     numXdim - integer - number of extra dimensions
!     numXvals - integer - size of valsXth
!     tunits - character*(*) - units for time, must be in the form of days
!      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
!      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
!     calendar - character*(*) - type of calendar.  Choose from 'noleap',
!      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
!      differ from modern.
!     * the following vars have an extra numXdim dimension (ignored if nameX='none') *
!     nameXth - character*(*) - name of Xth dimension
!     longXth - character*(*) - long name for Xth dimension variable
!     unitsXth - character*(*) - units for Xth dimension variable
!     axisXth - character*(*) - axis for Xth dimension variable
!     nXth - integer - size of Xth dimension
!     valsXth - real(nX) - values along Xth dimension
!     TODO: pos3rd variable is unused
!OUTPUT
!     idies - integer - id number of new file for later use
!     ierror - integer - error code, 0 = no error, < 0 = an error occured

      integer idies, nlon, nlat, ierror
      integer numXdim, numXvals, nXth(numXdim)
      character*(*) filen, title, source, history, tunits, calendar
      character*(*), dimension(numXdim) :: nameXth, longXth, unitsXth, axisXth
      real*8 alons(nlon), alats(nlat)
      real*8 valsXth(numXvals,numXdim)

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

      integer ierr, iddlon, idlon, iddlat, idlat
      integer, dimension(numXdim) :: iddXth, idXth
      integer iddtime, idtime, idtw, ll
      integer idims(2)

      ierror = 0

!     Open file
!     ---------
      ! if compression is requested, create a netcdf-4 file
      if ( env_compressout .eq. 0 ) then
         ierr = NF_CREATE(filen,NF_CLOBBER,idies)
      else
         ierr = NF_CREATE(filen,OR(NF_CLOBBER, OR(NF_NETCDF4, NF_CLASSIC_MODEL)),idies)
      end if
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 1 in inifile_dims, filen= '//filen
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

!     Define global attributes
!     ------------------------
      ll = lenchrib(title)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll+1,&
      title(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 2 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchrib(source)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll+1,&
      source(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 3 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchrib(history)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll+1,&
      history(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 4 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchrib(calendar)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll+1,&
      calendar(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 5 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'conventions',9,&
      'NCAR-CSM'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 6 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

!     Define dimensions
!     -----------------
      ierr = NF_DEF_DIM(idies,'longitude',nlon,iddlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 7 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'latitude',nlat,iddlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 8 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      do i = 1, numXdim
         iddXth(i) = -999
         if (nameXth(i) .ne. 'none') then
            ierr = NF_DEF_DIM(idies,nameXth(i),nXth(i),iddXth(i))
            if (ierr .ne. NF_NOERR) then
               print *, 'Error 9 in inifile_dims'
               print *, NF_STRERROR(ierr)
               ierror = -1
               return
            end if
         end if
      end do ! numXdim
     ierr = NF_DEF_DIM(idies,'time',NF_UNLIMITED,iddtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 10 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

!     Define dimension variables
!     --------------------------
      ierr = NF_DEF_VAR(idies,'longitude',NF_DOUBLE,1,iddlon,idlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 11 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'latitude',NF_DOUBLE,1,iddlat,idlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 12 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      do i = 1, numXdim
         if (iddXth(i) .ne. -999) then
            ierr = NF_DEF_VAR(idies,nameXth(i),NF_DOUBLE,1,iddXth(i),idXth(i))
            if (ierr .ne. NF_NOERR) then
               print *, 'Error 13 in inifile_dims'
               print *, NF_STRERROR(ierr)
               ierror = -1
               return
            end if
         end if
      end do ! numXdim
      ierr = NF_DEF_VAR(idies,'time',NF_DOUBLE,1,iddtime,idtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 14 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

!     Define variables time_weights and date
!     --------------------------------------
!     removed these vars - useless

!     Attributes for dimension variables
!     ----------------------------------
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'long_name',10,&
      'longitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 15 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'units',13,&
      'degrees_east'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 16 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'long_name',9,&
      'latitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 17 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'units',14,&
      'degrees_north'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 18 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      do i = 1, numXdim
         if (iddXth(i) .ne. -999) then
            ll = lenchrib(longXth(i))
            ierr = NF_PUT_ATT_TEXT(idies,idXth(i),'long_name',&
                 ll+1,longXth(i)(1:ll)//char(0))
            if (ierr .ne. NF_NOERR) then
               print *, 'Error 19 in inifile_dims'
               print *, NF_STRERROR(ierr)
               ierror = -1
               return
            end if
            ll = lenchrib(unitsXth(i))
            ierr = NF_PUT_ATT_TEXT(idies,idXth(i),'units',&
                 ll+1,unitsXth(i)(1:ll)//char(0))
            if (ierr .ne. NF_NOERR) then
               print *, 'Error 20 in inifile_dims'
               print *, NF_STRERROR(ierr)
               ierror = -1
               return
            end if
            ll = lenchrib(axisXth(i))
            ierr = NF_PUT_ATT_TEXT(idies,idXth(i),'axis',&
                 ll+1,axisXth(i)(1:ll)//char(0))
            if (ierr .ne. NF_NOERR) then
               print *, 'Error 21 in inifile_dims'
               print *, NF_STRERROR(ierr)
               ierror = -1
               return
            end if
         end if
      end do ! numXdim
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'long_name',5,&
      'time'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 22 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchrib(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll+1,&
      tunits(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 23 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

!     End define mode, enter data mode
!     --------------------------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 24 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

!     Put dimension variables except for time
!     ---------------------------------------
      ierr = NF_PUT_VAR_DOUBLE(idies,idlon,alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 25 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_DOUBLE(idies,idlat,alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 26 in inifile_dims'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      do i = 1, numXdim
         if (iddXth(i) .ne. -999) then
            ierr = NF_PUT_VAR_DOUBLE(idies,idXth(i),valsXth(:,i))
            if (ierr .ne. NF_NOERR) then
               print *, 'Error 27 in inifile_dims'
               print *, NF_STRERROR(ierr)
               ierror = -1
               return
            end if
         end if
      end do ! numXdim

!     Don't close file
!     ----------------
      return
end subroutine inifile_dims

! original version which initializes a file with one extra dimension only
subroutine inifile(idies,filen,title,source,history,nlon,alons,nlat,alats, &
                   name3rd,long3rd,units3rd,axis3rd,n3rd,vals3rd,pos3rd, &
                   tunits,calendar,ierror)

      integer idies, nlon, nlat, n3rd, ierror
      character*(*) filen, title, source, history, name3rd
      character*(*) long3rd, units3rd, pos3rd, axis3rd, tunits, calendar
      real*8 alons(nlon), alats(nlat), vals3rd(n3rd)

      ! extra dim variables passed to inifile
      integer, parameter :: numXdim = 1
      integer nXth(numXdim)
      character*80, dimension(numXDim) :: nameXth, longXth, unitsXth, axisXth
      real*8 valsXth(n3rd,numXdim)

      ! fill variables passed to inifile_dims
      nameXth(1)   = name3rd
      longXth(1)   = long3rd
      unitsXth(1)  = units3rd
      axisXth(1)   = axis3rd
      nXth(1)      = n3rd
      valsXth(1:n3rd,1) = vals3rd(:)

      call inifile_dims(idies,filen,title,source,history,nlon,alons,nlat,alats, &
                        nameXth,longXth,unitsXth,axisXth,nXth,valsXth,1,n3rd, &
                        tunits,calendar,ierror)

      return

end subroutine inifile
