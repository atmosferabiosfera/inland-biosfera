#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! if --enable-highres_map was not passed to configure, do not compile
#ifdef HRMAP

! ---------------------------------------------------------------------
subroutine whrvegmap (nday)
! ---------------------------------------------------------------------
! writes out hrmap files
! ---------------------------------------------------------------
      use inland_parameters
      use inland_control, only: iyear, iyear0, outdir, env_debug
      use inland_combcs
      use inland_subgrid

      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

! -----------------------------------------------------------------------
! input variables
      integer :: nday      ! number of days run since iyear0

! local variables
      integer :: idies
      integer :: mstep, istat, i, j, k, zero
      integer :: istart(4), icount(4)  ! for writing vars
      real*8 :: pindex(npft)    ! index used for pfts and canopies
      character*8 tmpdate       ! temp. variable to hold date
      character*8 cdate         ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*21 tunits       ! time units
      character*80 dimnames(4)  ! names of dimensions for vars
      character*1024 filen      ! file name
      character*4 chyear
      real*8 ftime,       & ! real form of nday
             tweight        ! number of days in yearly average

      ! use a real*8 dummy array for code simplicity (no need for writevarint)
      ! but this uses up more memory (50%) than using an int array
      !integer, dimension(:), allocatable :: cdummy2
      real*8, dimension(:), allocatable :: cdummy3
      integer dummyval
      
      logical save_ihrtileparent

      data istart / 1,1,1,1 /,  &
           icount / 1,1,1,1 /
! ---------------------------------------------------------------------
! FIXME: no point in using 'data' above if we change it right away to a
!       variable that is not compatible with 'data' assignment.
      icount(1) = hrmapnlon
      icount(2) = hrmapnlat

      zero = 0 ! you kiddin!

!      allocate(cdummy2(hrmapnlon*hrmapnlat))
      allocate(cdummy3(hrmapnlon*hrmapnlat))

! current time value and step: make ftime Jan 1 of this year
! instead of Dec 31
      ftime = nday - ndaypy + 1
      tweight = ndaypy
!     mstep = iyear - iyear0 + 1
      mstep = 1

      if (iyear .lt. 10) then
         write(chyear(1:1),'(i1)') zero
         write(chyear(2:2),'(i1)') zero
         write(chyear(3:3),'(i1)') zero
         write(chyear(4:4),'(i1)') iyear
      elseif (iyear .lt. 100) then
         write(chyear(1:1),'(i1)') zero
         write(chyear(2:2),'(i1)') zero
         write(chyear(3:4),'(i2)') iyear
      elseif (iyear .lt. 1000) then
         write(chyear(1:1),'(i1)') zero
         write(chyear(2:4),'(i3)') iyear
      else
         write(chyear,'(i4)') iyear
      endif

! tdate is ANN (3 char) followed by this year (4 char)
      tdate='ANN0000'//char(0)//char(0)//char(0)
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      elseif (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      elseif (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if

! first time only
      if (mstep .eq. 1) then
         call date_and_time(tmpdate)
         cdate(1:2) = tmpdate(5:6)
         cdate(3:3) = '/'
         cdate(4:5) = tmpdate(7:8)
         cdate(6:6) = '/'
         cdate(7:8) = tmpdate(3:4)

! time units is days since Dec 31 of the year before iyear0
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1

! define plant functional types, canopies with indices
         do 100 i = 1, npft
            pindex(i) = i
100      continue


      end if !mstep .eq. 1

      filen = trim(outdir)//'/inland-hrmap-'//chyear//'.nc'

      if ( env_debug.gt.1 ) print *,'creating '//trim(filen)

      if (myid .eq. 0) then
         if (mstep .eq. 1) then
            call inifile(idies,filen,'Inland', &
                         'inland whrvegmap',cdate,hrmapnlon,hrmaplonvalues,&
                         hrmapnlat,hrmaplatvalues,'none','none','none','Z', &
                         1,pindex,'none',tunits,'gregorian',istat)

! inivar

            dimnames(1) = 'longitude'
            dimnames(2) = 'latitude'
            dimnames(3) = 'time'

            save_ihrtileparent = .FALSE.
            save_ihrtileparent = .TRUE.
            if ( save_ihrtileparent ) then
               call inivar_type(idies,'ihrtileparent','parent tile index','none',3,dimnames, &
                    real(NF_FILL_SHORT,8),istat,NF_SHORT)
            end if
            call inivar_byte(idies,'vegtype','parent tile vegtype','none',3,dimnames,istat)
            !call inivar_byte(idies,'landuse','parent tile landuse','none',3,dimnames,istat)

            call closefile(idies,istat)

          end if

         call openfile(idies,filen,istat)

      endif

! loop over ihrtileparent for each variable - costs more time but less memory

! use writevar_notile for now, because we don't need the complexity (tiles, copying, etc) of writevar
! eventually, only need ihrtileparent, because we can make any maps using this var and inland-hrmap

! ihrtileparent
      if ( save_ihrtileparent ) then
      cdummy3(:) = NF_FILL_SHORT
      k = 0
      do j = 1, hrmapnlat
         do i = 1, hrmapnlon
            k = k + 1
            ! ihrtileparent is integer*2 with a -32768 offset
            !cdummy3(k) = ihrtileparent(j,i) + ihrtileparentoffset
            dummyval = ihrtileparent(j,i) + ihrtileparentoffset
            if ( dummyval .ne. 0 ) cdummy3(k) = dummyval
         end do
      end do

      if (myid .eq. 0) then
         istart(3) = mstep
         icount(3) = 1
         call writevar_notile(filen,idies,'ihrtileparent',cdummy3,istart,icount,ftime,tweight, &
                       tdate,istat)
         if (istat .ne. 0) then
            write(*,*) 'ERROR in whrvegmap, ihrtileparent'
            goto 999
         end if
      end if
      end if

! vegtype
      cdummy3(:) = NF_FILL_BYTE
      k = 0
      do j = 1, hrmapnlat
         do i = 1, hrmapnlon
            k = k + 1
            dummyval = ihrtileparent(j,i) + ihrtileparentoffset
            if ( dummyval .ne. 0 ) cdummy3(k) = xinveg(dummyval)
         end do
      end do

      if (myid .eq. 0) then
         istart(3) = mstep
         icount(3) = 1
         call writevar_notile(filen,idies,'vegtype',cdummy3,istart,icount,ftime,tweight, &
                       tdate,istat)
         istat = 0
         if (istat .ne. 0) then
            write(*,*) 'ERROR in whrvegmap, vegtype'
            goto 999
         end if
      endif

! landuse
#if 0
      cdummy2(:) = NF_FILL_BYTE
      k = 0
      do j = 1, hrmapnlat
         do i = 1, hrmapnlon
            k = k + 1
            dummyval = ihrtileparent(j,i) + ihrtileparentoffset
            if ( dummyval .ne. 0 ) cdummy2(k) = landuse(dummyval)
         end do
      end do

      if (myid .eq. 0) then
         istart(3) = mstep
         icount(3) = 1
         call writevar_notile(filen,idies,'landuse',cdummy2,istart,icount,ftime,tweight, &
                       tdate,istat)
         istat = 0
         if (istat .ne. 0) then
            write(*,*) 'ERROR in whrvegmap, landuse'
            goto 999
         end if
      endif
#endif

!      deallocate(cdummy2)
      deallocate(cdummy3)

      call closefile(idies,istat)

   return
   
999   write(*,*) ''
   stop


end subroutine whrvegmap

#endif
