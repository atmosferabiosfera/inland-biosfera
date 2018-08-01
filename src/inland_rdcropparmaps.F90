#include "inland_config.h"


  
! ---------------------------------------------------------------------
subroutine rdcropparmaps(iwest,jnorth)
! ---------------------------------------------------------------------
! gabriel abrahao: Reads yearly NetCDF maps of crop-related parameters
!


      use inland_parameters
      use inland_combcs
      use inland_control, only: iyear, env_debug, datadir
      use inland_comwork
      use inland_comcrop

      implicit none

      integer i,      & ! loop indice on years after istyr
              istat,  & ! error flag for netcdf
              istyr,  & ! 1st year in data files
              ndim,   & ! number of dimensions
              ntime,  & ! loop indices
              iwest,  & ! 1st lon index for subset
              jnorth    ! 1st lat index for subset

      character*1024 filen, suffix
      character*1024 :: subdirectory='crop_params'
      integer istart(2), icount(2)
      
      logical :: file_e

      real*8,dimension(:),allocatable :: plantdoy

      data istart / 1,1 /, icount / nlon,nlat /

      allocate(plantdoy(lbeg:lend))
      plantdoy = 0.0d0

      !gabriel apagar
      write(*,*) "======================================================================================== PASSEI inland_rdcropparmaps.F90"
      write(*,*) "shape(pdmin) = ",shape(pdmin)
      write(*,*) "pmmin(1,:) = ",pmmin(1,:)
      write(*,*) "pmmin(3,:) = ",pmmin(3,:)
      write(*,*) "datadir = ",datadir

      write(*,*) "shape(plantdoy) = ",shape(plantdoy)
      write(*,*) "plantdoy(1) = ",plantdoy(1)
      write(*,*) "ndaypm(2) = ",ndaypm(2)
      write(*,*) "mon_from_doy(274.0d0) = ",mon_from_doy(274.0d0)
      write(*,*) "mon_from_doy(293.17d0) = ",mon_from_doy(293.17d0)
      ! write(*,*) "mon_from_doy(182) = ",mon_from_doy(182)


      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

      ! icount(3) = 1
      ! icount(4) = 1

! Planting DOY maps for each year, still assuming a single crop is planted

      write(suffix,'(A,I4,A)') 'plantdoy.',iyear,'.nc'
      filen = trim(datadir)//'/'//trim(subdirectory)//'/'//trim(suffix)

      if ( env_debug .gt. 0 ) print *,'reading planting doy from '//trim(filen)

      !gabriel apagar
      write(*,*)"trim(filen) = ",trim(filen)

      ! make sure this file exists, if not print error and exit
      inquire( file=trim(filen), exist=file_e )
      if ( .not. file_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input file '//trim(filen)//' does not exist!'
         stop 1
      end if

      aname = 'plantdoy'
      call readvar(filen,aname,plantdoy,istart,icount,0,istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading plantdoy'
         stop 1
      end if

      !gabriel apagar
      write(*,*) "shape(plantdoy) = ",shape(plantdoy)
      write(*,*) "plantdoy(1) = ",plantdoy(1)
      write(*,*) "plantdoy(3) = ",plantdoy(3)
      write(*,*) "icroptype = ",icroptype

      !Replace pdmin and pmmin from icroptype with the read data, converted to day and month
      do i = lbeg,lend
        pdmin(i,icroptype) = day_from_doy(plantdoy(i))
        pmmin(i,icroptype) = mon_from_doy(plantdoy(i))       
      end do

      !gabriel apagar
      write(*,*) "======================================================================================== LEU"
      write(*,*) "shape(pdmin) = ",shape(pdmin)
      write(*,*) "pmmin(1,:) = ",pmmin(1,:)
      write(*,*) "pmmin(3,:) = ",pmmin(3,:)


! return to main program
      return

9000  format (1x,'ERROR in subroutine rdcropparmaps')
9010  format (1x,' ')
9020  format (1x,'number of land points: ', i10)


contains ! Functions to convert doy to day and month
  real*8 function day_from_doy(doy)
    real*8 doy
    integer i,j,mon
    i = 1
    j = 1
    mon = 1
    do while (i.lt.doy)
      if (j.eq.ndaypm(mon)) then
        mon = mon + 1
        j = 0
      end if
      j = j+1
      i = i+1
    end do


    day_from_doy=j
    return
  end function day_from_doy

  real*8 function mon_from_doy(doy)
    real*8 doy
    integer i,j,mon
    i = 1
    j = 1
    mon = 1
    do while (i.lt.doy)
      if (j.eq.ndaypm(mon)) then
        mon = mon + 1
        j = 0
      end if
      j = j+1
      i = i+1
    end do


    mon_from_doy=mon
    return
  end function mon_from_doy

end subroutine rdcropparmaps

