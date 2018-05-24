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
! These subroutines can be used as a primary interface to the ies
! format netcdf files.  You can also use the slightly lower level functions
! found in ies.f or the low level netcdf commands (see the netcdf V3
! users guide).
! SUBROUTINE LIST:
! inifile - create a new file with dimensions and variables for
! longitude, latitude, an optional 3rd real dimension, and time
! Time must have units of days since a date!!!
! inifilec - create a new file with dimensions and variables for
! longitude, latitude, a 3rd character-variable dimension, and time
! Time must have units of days since a date!!!
! inivar - initialize a new variable
! endini - end initialization of file: sync file to disk, close file
! writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
! readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
! .....................................................................
! inivar - initialize a new variable
! .....................................................................

! subroutine which initializes a variable with given type and missing value
subroutine inivar_type(idies,varname,longname,units,ndims,dimnames,valmissing, &
                  ierror,xtype)
!
      use inland_control, only: env_compressout, env_chunkout, env_floatout
!
      implicit none
!
!
! INPUT
! idies - integer - id number of a new file from a previous call
! to inifile, inifilec, or iesopen
! varname - character*(*) - name for new variable
! longname - character*(*) - descriptive name for variable
! units - character*(*) - units for variable, compatible with udunits
! ndims - integer - number of dimensions for this variable
! dimnames - character*(*)(ndims) - name of each dimension, in order
! valmissing - real - missing value, ignored if valmissing=0.
! xtype - netcdf data type to use (NF_INT, NF_DOUBLE, etc.)
! OUTPUT
!     ierror - integer - error code, 0 = no error, < 0 = error occurred

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

      integer :: idies, ndims, ierror, xtype
      character*(*) varname, longname, units, dimnames(ndims)
      real*8 :: valmissing
      real*4 :: valmissing2
      integer chunksize(ndims)
      integer :: ierr, idvar, iddims(5), i, lenchrib

      ierror = 0

! Put into redefine mode, but don't fail if already in define mode
! ----------------------------------------------------------------
      ierr = NF_REDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in inivar_type ',varname,idies
         print *, NF_STRERROR(ierr)
      end if

! Find id's of dimensions
! -----------------------
      do 100 i = 1, ndims
         ierr = NF_INQ_DIMID(idies,dimnames(i),iddims(i))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error 1 in inivar_type ',varname
            print *, i,idies,dimnames(i),iddims(i)
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
100   continue

! Define variable
! ---------------
      ierr = NF_DEF_VAR(idies,varname,xtype,ndims,iddims,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 2 in inivar_type ',varname,valmissing,xtype
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

! Enable DEFLATE compression if requested
! don't stop on error as it's not fatal
! chunking by lines is necessary due to gdal bug...
! ---------------
#ifdef NETCDF4
      if ( env_compressout .gt. 0 ) then
         ierr = NF_DEF_VAR_DEFLATE(idies,idvar,0,1,env_compressout)
         if (ierr .ne. NF_NOERR) then
            print *, 'Warning in inivar_type at NF_DEF_VAR_DEFLATE ',varname
            print *, NF_STRERROR(ierr)
         else if ( env_chunkout .eq. 1 ) then
            ! set chunksize to 1, except for longitude
            chunksize(:) = 1
            ! get lon length, assuming lon is first dim
            ierr = NF_INQ_DIMLEN(idies,iddims(1),i)
            if (ierr .ne. NF_NOERR) then
               print *, 'Warning in inivar_type at NF_INQ_DIMLEN ',varname
               print *, NF_STRERROR(ierr)
            else
               chunksize(iddims(1)) = i
               ierr = NF_DEF_VAR_CHUNKING(idies,idvar,NF_CHUNKED,chunksize)
               if (ierr .ne. NF_NOERR) then
                  print *, 'Warning in inivar_type at NF_DEF_VAR_CHUNKSIZE ',varname
                  print *, NF_STRERROR(ierr)
               end if              
            end if
         end if         
      end if
#else
      if ( env_compressout .gt. 0 ) then
         print *, 'Warning in inivar_type at NF_DEF_VAR_DEFLATE ',varname
         print *, 'netcdf-4 not available so you cannot use compression'
         env_compressout = 0
         env_chunkout = 0
      end if
#endif

! Define attributes
! -----------------
      ierr = NF_PUT_ATT_TEXT(idies,idvar,'long_name',lenchrib(longname), &
                             longname)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 3 in inivar_type ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idvar,'units',lenchrib(units),units)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 4 in inivar_type ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (valmissing .ne. 0.) then
         ! add missing value of same type as variable - watch out for overflows!
         if ( xtype .eq. NF_INT ) then
            ierr = NF_PUT_ATT_INT(idies,idvar,'missing_value',NF_INT,1,nint(valmissing))
         else if ( xtype .eq. NF_SHORT ) then
            ierr = NF_PUT_ATT_INT(idies,idvar,'missing_value',NF_SHORT,1,nint(valmissing))
         else if ( xtype .eq. NF_BYTE ) then
            ierr = NF_PUT_ATT_INT(idies,idvar,'missing_value',NF_BYTE,1,nint(valmissing))
         ! if user set INLAND_FLOATOUT, use NF_FLOAT type instead of NF_DOUBLE
         else if ( env_floatout .eq. 1 ) then 
            valmissing2 = valmissing !convert missing value to real*4
            ierr = NF_PUT_ATT_REAL(idies,idvar,'missing_value',NF_FLOAT,1,valmissing2)
         ! default to double with NF_FILL_DOUBLE missing value        
         else 
            ierr = NF_PUT_ATT_DOUBLE(idies,idvar,'missing_value',NF_DOUBLE,1,valmissing)
         end if
         if (ierr .ne. NF_NOERR) then
            print *, 'Error 5 in inivar_type ',varname
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

! Exit define mode
! ----------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error 6 in inivar_type ',varname
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      return
end subroutine inivar_type

! initializes a double (real*8) variable
! original version which calls inivar_type() with xtype=NF_DOUBLE/NF_FLOAT
! TODO remove nodata
subroutine inivar(idies,varname,longname,units,ndims,dimnames, &
                  ierror)

      use inland_control, only: env_floatout
      use inland_comwork, only: ocean

      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

      integer :: idies, ndims, ierror
      character*(*) varname, longname, units, dimnames(ndims)

      if ( env_floatout .eq. 1 ) then ! if user set INLAND_FLOATOUT, use NF_FLOAT
         call inivar_type(idies,varname,longname,units,ndims,dimnames,ocean, &
              ierror,NF_FLOAT)
      else ! default to double
         call inivar_type(idies,varname,longname,units,ndims,dimnames,ocean, &
              ierror,NF_DOUBLE)
      end if
         
end subroutine inivar

! initializes a byte (integer*1) variable
! new version which calls inivar_type() with xtype=NF_BYTE
subroutine inivar_byte(idies,varname,longname,units,ndims,dimnames,ierror)

      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

      integer :: idies, ndims, ierror, valmissing
      character*(*) varname, longname, units, dimnames(ndims)

      call inivar_type(idies,varname,longname,units,ndims,dimnames,&
                  real(NF_FILL_BYTE,8),ierror,NF_BYTE)
end subroutine inivar_byte

! initializes a integer variable
! new version which calls inivar_type() with xtype=NF_INT
subroutine inivar_int(idies,varname,longname,units,ndims,dimnames,ierror)

      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

      integer :: idies, ndims, ierror, valmissing
      character*(*) varname, longname, units, dimnames(ndims)

      call inivar_type(idies,varname,longname,units,ndims,dimnames,&
                  real(NF_FILL_INT,8),ierror,NF_INT)
end subroutine inivar_int

