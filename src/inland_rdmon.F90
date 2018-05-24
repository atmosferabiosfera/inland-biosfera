#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine rdmon(istyr, iwest, jnorth)

      use inland_parameters
      use inland_control, only: imonth, iyear, datadir
      use inland_combcs
      use inland_comwork

      implicit none

      integer istat, &  ! error flag for netcdf
              istyr,  & ! 1st year in data files - unused
              ndim, & !number of dimensions
              ntime,& ! loop indices
              iwest,  & ! 1st lon index for subset
              jnorth    ! 1st lat index for subset! number of dimensions

      character*1024 filen
      character*1024 :: directory
      character*80 :: suffix
      integer istart(4), icount(4)
      logical :: file_e

      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

      icount(3) = 1
      icount(4) = 12

      directory = trim(datadir)//'/monthly/'

      aname = 'wetd'
      write(suffix,'(A,I4,A)') '.monthly.',iyear,'.nc'
      filen = trim(directory)//trim(aname)//trim(suffix)

      ! make sure this file exists, if not print error and exit
      inquire( file=trim(filen), exist=file_e )
      if ( .not. file_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input file '//trim(filen)//' does not exist!'
         write (*,*) 'make sure INLAND_DATADIR is set to proper path and that file exists'
         stop 1
      end if

      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999

      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obswet(1,ntime))
      end do

      aname = 'temp'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obst(1,ntime))
      end do

      aname = 'trange'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obstrng(1,ntime))
      end do

      aname = 'prec'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obsprec(1,ntime))
      end do

      aname = 'wspd'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obswind(1,ntime))
      end do

      aname = 'cld'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obscld(1,ntime))
      end do

      aname = 'rh'
      filen = trim(directory)//trim(aname)//trim(suffix)
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),obsq(1,ntime))
      end do

! copy all 5 climatology fields to clm+anom fields for spin up
      ndim = npoi*12
      call scopya(ndim, obst, xintmon)
      call scopya(ndim, obswind, xinwindmon)
      call scopya(ndim, obstrng, xintrngmon)
      call scopya(ndim, obsprec, xinprecmon)
      call scopya(ndim, obscld, xincldmon)
      call scopya(ndim, obsq, xinqmon)
      call scopya(ndim, obswet, xinwetmon)

! return to main program
      return

9000  format (1x,'ERROR in subroutine rdmon')

9999  write(*,9000)
      print *, 'while reading '//trim(filen)
      stop 1

end subroutine rdmon
