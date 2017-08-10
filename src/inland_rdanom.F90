#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D IBIS model option."
#endif

! ---------------------------------------------------------------------
subroutine rdanom(iyranom,nanom,iy2,istyr,iwest,jnorth)
! ---------------------------------------------------------------------
!
! reads in anomalies for imonth+1
!
! this subroutine reads the monthly anomaly values of:
!
! temp - mean temperature
! trng - temperature range (not read in yet, saved for future use)
! prec - precipitation rate
! cld  - cloudiness
! rh   - relative humidity
! wspd - wind speed (not read in yet, saved for future use)
! wetd - wet days per month (not read in yet, saved for future use)
!
! and adds them to the climatological value for the month
      use inland_parameters
      use inland_control, only: imonth, iyear
      use inland_combcs
      use inland_comveg
      use inland_comwork

      implicit none

! Arguments
      integer iyranom,& ! year to start reading anoms
              iy2,    & ! final year of the run
              istyr,  & ! 1st year in data files
              iwest,  & ! 1st lon index for subset
              jnorth, & ! 1st lat index for subset
              nanom     ! # of years in the anomaly files

! Local variables
      real*8, dimension(:), pointer :: anom

      integer imon,  & ! month + 1
              iyr,   & ! number of years after begin of data file(istyr)
              istat, & ! error flag for netcdf
              i,     & ! loop indice on land points
              jyear    ! iyr divided by nanom (for looping through anomaly files)

      integer istart(4), icount(4) ! for reading rdanom vars
      character*80 filen

! ---------------------------------------------------------------------
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

! Point 'anom' to cdummy structure
      anom => cdummy

! determine which (if any) month to read
! If prior to November of the year before an anomaly year, then return.
! Else read anomalies for the following month.
      if (iyear .lt. iyranom-1) then
         return
      else if (iyear .eq. iyranom-1 .and. imonth .lt. 11) then
         return

!    If timestep equals December of final year in anomaly file and also equals
! the final year of the run, then set January anomalies to zero and return.
! If not the final year of the run, then loop back to start of anomaly file
! (see jyear).
      else if (iyear .eq. istyr+nanom-1 .and. imonth .eq. 12) then
         if (iyear .eq. iy2) then
            print *, 'WARNING: last month of run; no anomalies for January'
            print *, 'Using climatologies for month year ='
            print *, imonth+1,iyear+1
            do 4 i = lbeg, lend
               xint(i,1) = clmt(i,1)
!             xintrng(i,1) = clmtrng(i,1)
               xinprec(i,1) = clmprec(i,1)
               xincld(i,1) = clmcld(i,1)
               xinq(i,1) = clmq(i,1)
!             xinwind(i,1) = clmw(i,1)
!             xinwet(i,1) = clmwet(i,1)
4           continue
            return
         end if
      end if

      iyr = iyear-istyr
      imon = imonth + 1
      if (imon .eq. 13) then
         imon = 1
         iyr = iyr + 1
      end if
      jyear = iyr/nanom
      istart(4) = (iyr - nanom*jyear)*12 + imon

      if (iyr.gt.0 .and. (iyr - nanom*jyear).eq.0) then
         print *, 'WARNING: Attempted to read past last month in anomaly file'
         print *, 'Looping back to the beginning of the file'
      end if

      if (istart(4) .gt. 0) then
         print *, 'rdanom reading month year step ='
         print *, imon,iyr+istyr-nanom*jyear,istart(4)
      else
         print *, 'WARNING, anomalies begin in year ',istyr
         print *, 'Not reading in anomalies for month year ='
         print *, imon,iyr+istyr
         return
      end if

!    dummy variable example, 4-d, whose 3rd dim (level) = 1
!
!     aname = 'dummyv'
!     filen = 'input/anom/dummyv.danom.nc'
!     call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
!                 cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
!     if (istat .ne. 0) then
!        write(*,*) 'ERROR in rdanom, dummyv'
!        stop 1
!     end if
!     call arr2vec( work, anom )
!     do 10 i = 1, npoi
!        xindummyv(i,imon) = max (clmtrng(i,imon) + anom(i), 0.1)
!10   continue
!
!    mean temperature
      aname = 'temp'
      filen = 'input/anom/temp.danom.nc'
      call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
                   cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, temp'
         stop 1
      end if
      call arr2vec( work, anom )
      do 5 i = lbeg, lend
         xint(i,imon) = clmt(i,imon) + anom(i)
5     continue
! 
!    temperature range
!
!     aname = 'trange'
!     filen = 'input/anom/trange.danom.nc'
!     call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
!                  cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
!     if (istat .ne. 0) then
!        write(*,*) 'ERROR in rdanom, trange'
!        stop 1
!     end if
!     call arr2vec( work, anom )
!     do 10 i = 1, npoi
!        xintrng(i,imon) = max (clmtrng(i,imon) + anom(i), 0.1)
!10   continue
!
!    precipitation rate
      aname = 'prec'
      filen = 'input/anom/prec.danom.nc'
      call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
                   cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, prec'
         stop 1
      end if
      call arr2vec( work, anom )
      do 15 i = lbeg, lend
         xinprec(i,imon) = clmprec(i,imon) + anom(i)
15    continue

!    cloudiness
      aname = 'cld'
      filen = 'input/anom/cld.danom.nc'
      call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
                   cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, cld'
         stop 1
      end if
      call arr2vec( work, anom )
      do 20 i = lbeg, lend
         xincld(i,imon) = clmcld(i,imon) + anom(i)
20    continue

!    relative humidity
      aname = 'rh'
      filen = 'input/anom/rh.danom.nc'
      call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
                   cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, rh'
         stop 1
      end if
      call arr2vec( work, anom )
      do 25 i = lbeg, lend
         xinq(i,imon) = clmq(i,imon) + anom(i)
25    continue

!    wind speed
!     aname = 'wspd'
!     filen = 'input/anom/wspd.danom.nc'
!     call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
!                  cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
!     if (istat .ne. 0) then
!        write(*,*) 'ERROR in rdanom, wspd'
!        stop 1
!     end if
!     call arr2vec( work, anom )
!     do 30 i = 1, npoi
!        xinwind(i,imon) = clmw(i,imon) + anom(i)
!30   continue
!
!    wet days
!
!     aname = 'wetd'
!     filen = 'input/anom/wetd.danom.nc'
!     call readvar(filen,aname,'level',istart,icount,work,cdummy(1),cdummy(nlonsub+1),&
!                  cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
!     if (istat .ne. 0) then
!        write(*,*) 'ERROR in rdanom, wetd'
!        stop 1
!     end if
!     call arr2vec( work, anom )
!     do 35 i = 1, npoi
!        xinwet(i,imon) = clmwet(i,imon) + anom(i)
!35   continue
      return
end subroutine rdanom
