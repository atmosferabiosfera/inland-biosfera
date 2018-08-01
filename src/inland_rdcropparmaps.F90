#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine rdcropparmaps(iwest,jnorth)
! ---------------------------------------------------------------------
! gabriel abrahao: Reads yearly NetCDF maps of crop-related parameters
!

      use inland_parameters
      use inland_combcs
      use inland_control, only: iyear, env_debug
      use inland_comveg,  only: disturbl
      use inland_comwork

      implicit none

      integer i,      & ! loop indice on years after istyr
              istat,  & ! error flag for netcdf
              istyr,  & ! 1st year in data files
              ndim,   & ! number of dimensions
              ntime,  & ! loop indices
              iwest,  & ! 1st lon index for subset
              jnorth    ! 1st lat index for subset

      character*1024 filen, suffix
      character*1024 :: directory='input/transition/'
      integer istart(4), icount(4)
      
      logical :: file_e

      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
     
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

      icount(3) = 1
      icount(4) = 1

! deforestation maps for each year

      write(suffix,'(A,I4,A)') 'Forest',iyear,'.nc'
      filen = trim(directory)//trim(suffix)

      if ( env_debug .gt. 0 ) print *,'reading transition from '//trim(filen)

      ! make sure this file exists, if not print error and exit
      inquire( file=trim(filen), exist=file_e )
      if ( .not. file_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input file '//trim(filen)//' does not exist!'
         stop 1
      end if

      aname = 'transition'
      call readvar(filen,aname,disturbl,istart,icount,0,istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading disturbl'
         stop 1
      end if

! return to main program
      return

9000  format (1x,'ERROR in subroutine rdcropparmaps')
9010  format (1x,' ')
9020  format (1x,'number of land points: ', i10)

end subroutine rdcropparmaps

