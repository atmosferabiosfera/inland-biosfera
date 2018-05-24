#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! This subroutine reads in daily fields..
! ---------------------------------------------------------------------
subroutine rdday(jday, istyr, iwest, jnorth)
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control, only: imonth, iyear, datadir
      use inland_combcs
      use inland_comwork

      implicit none

! Arguments
      integer jday,   & ! day of the year
              istyr,  & ! 1st year in data files
              iwest,  & ! 1st lon index for subset
              jnorth    ! 1st lat index for subset

! Local variables
      integer i,    & ! loop indice on years after istyr
              istat   ! error flag for netcdf
      integer dimsiz  ! Gabriel Abrahao: number of dimensions in input file
      integer cnt     ! Gabriel Abrahao: loop index

      character*1024 filen
      character*1024 :: directory
      character*80 :: suffix
      integer, dimension(:), allocatable :: istart, icount
      logical :: file_e

     ! Gabriel Abrahao: in original 3d file was:
     ! integer istart(3), icount(3)
     ! data istart / 1,1,1 /, icount / nlon,nlat,1 /
      
      if ( rddaydims .eq. 0 ) then
        dimsiz = 4
      else
        dimsiz = 3
      end if          

      allocate(istart(dimsiz),icount(dimsiz))

      do cnt = 1,dimsiz
        istart(cnt) = 1        
      enddo

      icount(1)=nlon ! Gabriel Abrahao: these two lines appear useless, as they are overwritten right below, keeping just for sure
      icount(2)=nlat
      icount(3)=1
      if ( rddaydims .eq. 0 ) then
        icount(4) = 1
      end if      
      
      

      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

      if (iyear .lt. istyr) then
         print *, 'daily data begins in year ', istyr
         print *, 'not reading in daily data'
         stop 1
      end if

! count how many days since beginning of daily data
! daily data begin on Jan 1, istyr
      if (iyear .eq. istyr) then
         istart(dimsiz) = jday
      else
         istart(dimsiz) = 0
         do 10 i = istyr, iyear-1
            istart(dimsiz) = istart(dimsiz) + 365
            if (mod(i,dimsiz) .eq.0) then
               if (mod(i,100).ne.0) then
                  istart(dimsiz) = istart(dimsiz) + 1
               else if (mod(i/100,dimsiz).eq.0) then
                  istart(dimsiz) = istart(dimsiz) + 1
               end if
            end if
10       continue
        istart(dimsiz) = jday
      end if

      directory = trim(datadir)//'/daily/'

! read daily precip
      aname = 'prec'
      write(suffix,'(A,I4,A)') '.daily.',iyear,'.nc'
      filen = trim(directory)//trim(aname)//trim(suffix)

      ! make sure this file exists, if not print error and exit
      inquire( file=trim(filen), exist=file_e )
      if ( .not. file_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input file '//trim(filen)//' does not exist!'
         write (*,*) 'make sure INLAND_DATADIR is set to proper path and that file exists'
         stop 1
      end if

      call readvar(filen,aname,xinprecd,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! read daily temp
      aname = 'temp'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,xintd,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! read daily trange
!      aname = 'trange'
!      filen = 'input/daily/trange.daily.nc'
!      call arr2vec( work, xintrngd(1) )

! read daily cloudiness
      aname = 'cld'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,xincldd,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! read daily windspeed
      aname = 'wspd'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,xinwindd,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! read daily srelative humidity
      aname = 'rh'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,xinqd,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! read tmax temperature
      aname = 'tmax'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,xintmaxd,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! read tmin temperature
      aname = 'tmin'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,xintmind,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

      return

9999  if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, '//trim(aname)
         ! TODO is this test really needed?
         if (iyear .gt. 1997) then
            print *, 'Attempted to read past last day in file?'
         end if
         stop 1
      end if

end subroutine rdday
