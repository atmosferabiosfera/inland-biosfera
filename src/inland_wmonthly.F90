#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
subroutine wmonthly(nday)
! ---------------------------------------------------------------------
! writes out monthly files
!---------------------------------------------------------------
      use inland_parameters
      use inland_control, only: imonth, iyear, iyear0, outdir, hrmapfile, env_debug, isimfire, isimveg
      use inland_comfire
      use inland_comsoi
      use inland_comwork
      use inland_comsum
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

!-----------------------------------------------------------------------
!     local variables
      integer mstep,  & ! time step in netcdf file
              n, k, idies, idiespft, istat, i, zero

      integer istart(5), icount(5) ! for writing vars
      real*8 :: pindex(npft)    ! index used for pfts and canopies
      real*8 :: tindex(mlpt+1) ! index used for tiles (last one is average/water)
      character*8 tmpdate      ! temp. variable to hold date
      character*8 cdate        ! date to use in history attribute in files
      character*13 canopies(2)  ! canopy definitions
      character*1024 canopies2  ! canopy definitions
      character*21 tunits      ! units for time
      character*80 dimnames(5) ! names of dimensions for vars
      character*60 pftdef(npft) ! plant functional type defs
      character*1024 pftdef2
      character*60 vegtypedef(15) ! vegtype defs (15 currently)
      character*1024 vegtypedef2
      character*1024 tiledef
      character*1024 filen, filenpft     ! file name
      character*3  chmon(12)   ! month abbreviations
      character*4 chyear
      integer ilpt             ! mlpt index
      integer ndims            ! number of dims (including tile, not including pft)
      integer mlpt1            ! number of tiles in output file, including average if mlpt > 1
      character*1 endline   ! endline character

! extra dim variables passed to inifile_dims
      integer, parameter :: numXdim = 2
      integer, parameter :: numXvals = 100 ! pick a big number to make sure we have enough space
      integer nXth(numXdim)
      character*80, dimension(numXDim) :: nameXth, longXth, unitsXth, axisXth
      real*8 valsXth(numXVals,numXdim)
      real*8 ftime             ! real form of nday

      data chmon / 'JAN','FEB','MAR','APR','MAY','JUN', &
                   'JUL','AUG','SEP','OCT','NOV','DEC' /
      data istart / 1,1,1,1,1 /,    &
           icount / nlon,nlat,1,1,1 /

! FIXME: Better drop 'data' statement than make it and just change it!
      icount(1) = nlonsub
      icount(2) = nlatsub

      zero = 0 ! you kiddin!
      endline = char(10) ! this is not read properly by ncdump<4.3 when using netcdf-4 files

! ---------------------------------------------------------------------
!     current time value, step, time weight

      ftime = nday
!     mstep = 12*(iyear-iyear0) + imonth
      mstep = imonth
!     mstep = 1
!      write(*,*) 'wmonthly, test: mstep = 1'

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

! how many dims?
      if ( mlpt .gt. 1 ) then
         ndims = 4
         mlpt1 = mlpt + 1
      else
         ndims = 3
         mlpt1 = 1
      end if

! first time only
      if (mstep.eq.1) then

! initialize snow layer indices, pft names, etc
         call date_and_time(tmpdate)
         cdate(1:2) = tmpdate(5:6)
         cdate(3:3) = '/'
         cdate(4:5) = tmpdate(7:8)
         cdate(6:6) = '/'
         cdate(7:8) = tmpdate(3:4)

! time units is days since Dec 31 of the year before iyear0
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1

! define plant functional types
         do i = 1, npft
            pindex(i) = i
         end do

! define tile dimensions
         do ilpt = 1, mlpt1
            tindex(ilpt) = ilpt
         end do
         ilpt = 1

      end if

! pft defs
         pftdef(1)  = ' 1: trbrevtr - tropical broadleaf evergreen trees'
         pftdef(2)  = ' 2: trbrdetr - tropical broadleaf drought-deciduous trees'
         pftdef(3)  = ' 3: wtbrevtr - warm-temperate broadleaf evergreen trees'
         pftdef(4)  = ' 4: tecoevtr - temperate conifer evergreen trees'
         pftdef(5)  = ' 5: tebrdetr - temperate broadleaf cold-deciduous trees'
         pftdef(6)  = ' 6: bocoevtr - boreal conifer evergreen trees'
         pftdef(7)  = ' 7: bocodetr - boreal conifer cold-deciduous trees'
         pftdef(8)  = ' 8: bobrdetr - boreal broadleaf cold-deciduous trees'
         pftdef(9)  = ' 9: evsh - evergreen shrubs'
         pftdef(10) = '10: desh - deciduous shrubs'
         pftdef(11) = '11: c4gr - warm (c4) grasses'
         pftdef(12) = '12: c3gr - cool (c3) grasses'
         write(pftdef2,'(A,A,A,A,A,A,A,A,A,A,A,A,A)') endline, &
              pftdef(1)//endline,pftdef(2)//endline,pftdef(3)//endline, &
              pftdef(4)//endline,pftdef(5)//endline,pftdef(6)//endline, &
              pftdef(7)//endline,pftdef(8)//endline,pftdef(9)//endline, &
              pftdef(10)//endline,pftdef(11)//endline,pftdef(12)//char(0)
         
! vegtype defs from inland_vegmap.F90
         vegtypedef(1)  = ' 1: tropical evergreen forest / woodland'
         vegtypedef(2)  = ' 2: tropical deciduous forest / woodland'
         vegtypedef(3)  = ' 3: temperate evergreen broadleaf forest / woodland'
         vegtypedef(4)  = ' 4: temperate evergreen conifer forest / woodland'
         vegtypedef(5)  = ' 5: temperate deciduous forest / woodland'
         vegtypedef(6)  = ' 6: boreal evergreen forest / woodland'
         vegtypedef(7)  = ' 7: boreal deciduous forest / woodland'
         vegtypedef(8)  = ' 8: mixed forest / woodland'
         vegtypedef(9)  = ' 9: savanna'
         vegtypedef(10) = '10: grassland / steppe'
         vegtypedef(11) = '11: dense shrubland'
         vegtypedef(12) = '12: open shrubland'
         vegtypedef(13) = '13: tundra'
         vegtypedef(14) = '14: desert '
         vegtypedef(15) = '15: polar desert / rock / ice'
         write(vegtypedef2,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)') endline, &
              vegtypedef(1)//endline,vegtypedef(2)//endline, &
              vegtypedef(3)//endline,vegtypedef(4)//endline, &
              vegtypedef(5)//endline,vegtypedef(6)//endline, &
              vegtypedef(7)//endline,vegtypedef(8)//endline, &
              vegtypedef(9)//endline,vegtypedef(10)//endline, &
              vegtypedef(11)//endline,vegtypedef(12)//endline, &
              vegtypedef(13)//endline,vegtypedef(14)//endline, &
              vegtypedef(15)//char(0)
         
! canopy defs
         canopies(1) = 'lower canopy'
         canopies(2) = 'upper canopy'
         canopies2 = '1: '//canopies(1)//' /  2: '//canopies(2)

! tiles defs
         write(tiledef, '(I1,A,I1,A,I1,A)') & 
              1, ' - ', mlpt, ' : tiles  /  ', mlpt+1, ' : area average'

! normalize the albedo calculation (not at each time step)
      do i = lbeg, lend
        amreflect(i) = amreflect(i) / (amsolar(i) + epsilon)
      end do

!
! Recording variables in a single file.  
!
      filen = trim(outdir)//'/inland-monthly-'//chyear//'.nc'
      filenpft = trim(outdir)//'/inland-monthly-pft-'//chyear//'.nc'

      if (env_debug .gt. 1 ) print *,'creating '//trim(filen)//' and '//trim(filenpft)

      if (myid .eq. 0) then
         if (mstep .eq. 1) then

               ! fill variables passed to inifile_dims
               nameXth(1)  = 'pft'
               longXth(1)  = 'plant functional type'
               unitsXth(1) = 'none'
               axisXth(1)  = 'E'
               nXth(1)     = npft
               valsXth(1:npft,1) = pindex
               nameXth(2)  = 'tile'
               longXth(2)  = 'tile'
               unitsXth(2) = 'none'
               axisXth(2)  = 'Z'
               nXth(2)     = mlpt1
               valsXth(:,2) = 0
               if ( mlpt1 .gt. 1 ) then
                  valsXth(1:mlpt1,2) = tindex
               else
                  valsXth(1:2,2) = tindex
               end if

               call inifile_dims(idies,filen,'Inland', &
                         'inland wmonthly',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                         tunits,'gregorian',istat)
               call inifile_dims(idiespft,filenpft,'Inland', &
                         'inland wmonthly pft',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                         tunits,'gregorian',istat)

! add global attribute to define pfts, vegtypes and canopies with text
 ! use netcdf low-level commands (must go back into redef mode)

            istat = NF_REDEF(idies)
            if ( mlpt .gt. 1 ) then
               istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'tile_def', &
                    LEN_TRIM(tiledef), trim(tiledef))
            end if
            istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def', &
                 LEN_TRIM(canopies2), canopies2)
            istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_def', &
                 LEN_TRIM(pftdef2), pftdef2)
            istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'vegtype_def', &
                 LEN_TRIM(vegtypedef2), vegtypedef2)
            istat = NF_ENDDEF(idies)
            if (istat .ne. NF_NOERR) then
               print *, 'Error in wyearly'
               print *, NF_STRERROR(istat)
            end if
            istat = NF_REDEF(idiespft)
            if ( mlpt .gt. 1 ) then
               istat = NF_PUT_ATT_TEXT(idiespft,NF_GLOBAL,'tile_def', &
                    LEN_TRIM(tiledef), trim(tiledef))
            end if
            istat = NF_PUT_ATT_TEXT(idiespft,NF_GLOBAL,'canopy_def', &
                 LEN_TRIM(canopies2), canopies2)
            istat = NF_PUT_ATT_TEXT(idiespft,NF_GLOBAL,'pft_def', &
                 LEN_TRIM(pftdef2), pftdef2)
            istat = NF_PUT_ATT_TEXT(idiespft,NF_GLOBAL,'vegtype_def', &
                 LEN_TRIM(vegtypedef2), vegtypedef2)
            istat = NF_ENDDEF(idiespft)
            if (istat .ne. NF_NOERR) then
               print *, 'Error in wyearly'
               print *, NF_STRERROR(istat)
            end if

! inivar

            dimnames(1) = 'longitude'
            dimnames(2) = 'latitude'
            dimnames(3:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(3) = 'tile'

            call inivar(idies,'aet','average evapotranspiration','mm/day',ndims, &
                        dimnames,istat)
            call inivar(idies,'transu','monthly upper canopy transpiration', &
                        'mm/day',ndims,dimnames,istat)
            call inivar(idies,'transl','monthly lower canopy transpiration', &
                        'mm/day',ndims,dimnames,istat)
            call inivar(idies,'suvap','monthly soil evaporation','mm/day',ndims, &
                        dimnames,istat)
            call inivar(idies,'invap','monthly interception loss','mm/day',ndims,&
                        dimnames,istat)
            call inivar(idies,'trunoff','average total runoff','mm/day',ndims,   &
                        dimnames,istat)
            call inivar(idies,'srunoff','average surface runoff','mm/day',ndims, &
                        dimnames,istat)
            call inivar(idies,'drainage','average drainage','mm/day',ndims, &
                        dimnames,istat)
            call inivar(idies,'vwc','average volumetric water content',  &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'awc','average plant available water content', &
                        'cm',ndims,dimnames,istat)
            call inivar(idies,'snod','average snow depth','meters',ndims,dimnames,&
                        istat)
            call inivar(idies,'snof','average snow fraction','m^2/m^3',ndims,    &
                        dimnames,istat)
            call inivar(idies,'albedo','average albedo',&
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'npptot','total npp of carbon', &
                        'kg m-2 month-1',ndims,dimnames,istat)
            call inivar(idies,'laiu','average lai upper canopy', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'lail','average lai lower canopy', &
                        'fraction',ndims,dimnames,istat)
            if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
               call inivar(idies,'Pbio','Fire prob due to biomass', &
                           'fraction',ndims,dimnames,istat)
               call inivar(idies,'Pmoi','Fire prob due to wetness', &
                           'fraction',ndims,dimnames,istat)
               call inivar(idies,'Pign','Fire prob due to ignition', &
                           'fraction',ndims,dimnames,istat)
               call inivar(idies,'Pfire','Fire prob', &
                           'fraction',ndims,dimnames,istat)
               call inivar(idies,'srate','fire spread rate','km/h',ndims,&
                           dimnames,istat)
               call inivar(idies,'abday','area burned in 1 day','km^2',ndims,&
                           dimnames,istat)
               call inivar(idies,'burnfrac','burned fraction','fraction',ndims,&
                           dimnames,istat)
               call inivar(idies,'burnarea','burned area','km^2',ndims,&
                           dimnames,istat)
               end if
            call inivar(idies,'transET','monthly average transpiration:ET ratio', &
                         'fraction',ndims,dimnames,istat)
            call inivar(idies,'amtotnleach','average total n leaching', &
                         'kg/ha/day',ndims,dimnames,istat)
            call inivar(idies,'no3leach','average nitrate leaching', &
                         '(kg/ha/day)',ndims,dimnames,istat)
            call inivar(idies,'tsoi','average soil temperature', &
                         'degC',ndims,dimnames,istat)
            call inivar(idies,'wisoi','fraction of soil pore space containing ice', &
                         'fraction',ndims,dimnames,istat)
            call inivar(idies,'wsoi','fraction of soil pore space containing liquid water', &
                         'fraction',ndims,dimnames,istat)
            call inivar(idies,'sens','average sensible heat flux', &
                         'W/m^2',ndims,dimnames,istat)
            call inivar(idies,'latent','average latent heat flux', &
                         'W/m^2',ndims,dimnames,istat)
! variables which have pft dimension
! for reading these files with cdo, when tiling is used, see "Subgrid Tiling" in docs/README.advanced

            dimnames(3) = 'pft'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'
            call inivar(idiespft,'burnpft','area burned in each pft','km^2',&
                        ndims+1,dimnames,istat)

            call closefile(idies,istat)
            call closefile(idiespft,istat)

         end if
         
         call openfile(idies,filen,istat)
         call openfile(idiespft,filenpft,istat)


      endif

      if (env_debug .gt. 1 ) print *,'writing to '//trim(filen)

      if ( mlpt .gt. 1 )  then
         istart(3) = 1
         icount(3) = mlpt1
      end if
      istart(ndims) = mstep
      icount(ndims) = 1

      call writevar(filen,idies,'aet',amaet,istart,icount,ftime,istat)
      call writevar(filen,idies,'transu',amtransu,istart,icount,ftime,istat)
      call writevar(filen,idies,'transl',amtransl,istart,icount,ftime,istat)
      call writevar(filen,idies,'suvap',amsuvap,istart,icount,ftime,istat)
      call writevar(filen,idies,'invap',aminvap,istart,icount,ftime,istat)
! trunoff, srunoff, drainage
      call writevar(filen,idies,'trunoff',amtrunoff,istart,icount,ftime,istat)
      call writevar(filen,idies,'srunoff',amsrunoff,istart,icount,ftime,istat)
      call writevar(filen,idies,'drainage',amdrainage,istart,icount,ftime,istat)
! soil temperature
! soil moisture, ice, volumetric water content, plant available water
      call writevar(filen,idies,'vwc',amvwc,istart,icount,ftime,istat)
      call writevar(filen,idies,'awc',amawc,istart,icount,ftime,istat)
! snow depth
      call writevar(filen,idies,'snod',amsnod,istart,icount,ftime,istat)
! snow fraction
      call writevar(filen,idies,'snof',amsnof,istart,icount,ftime,istat)
! solar radiation
! albedo
      call writevar(filen,idies,'albedo',amreflect,istart,icount,ftime,istat)
! downward and upward infrared radiation
! sensible heat flux
! latent heat flux
! total net primary productivity
      call writevar(filen,idies,'npptot',amnpptot,istart,icount,ftime,istat)
! leaf area index upper and lower
      call writevar(filen,idies,'laiu',amlaiu,istart,icount,ftime,istat)
      call writevar(filen,idies,'lail',amlail,istart,icount,ftime,istat)
! co2ratio
      call writevar(filen,idies,'transET',amtratio,istart,icount,ftime,istat)
      call writevar(filen,idies,'amtotnleach',amtotnleach,istart,icount,ftime,istat)
      call writevar(filen,idies,'no3leach',amno3leach,istart,icount,ftime,istat)
      call writevar(filen,idies,'tsoi',amtsoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'wisoi',wisoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'wsoi',amwsoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'sens',amsens,istart,icount,ftime,istat)
      call writevar(filen,idies,'latent',amlatent,istart,icount,ftime,istat)

! fire vars
      if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
         call writevar(filen,idies,'Pbio',ampbio,istart,icount,ftime,istat)
         call writevar(filen,idies,'Pmoi',ampmoi,istart,icount,ftime,istat)
         call writevar(filen,idies,'Pign',ampign,istart,icount,ftime,istat)
         call writevar(filen,idies,'Pfire',ampfire,istart,icount,ftime,istat)
         call writevar(filen,idies,'srate',amsrate,istart,icount,ftime,istat)
         call writevar(filen,idies,'abday',amabday,istart,icount,ftime,istat)
         call writevar(filen,idies,'burnfrac',amburnfrac,istart,icount,ftime,istat)
         call subgrid_calculate_area(amburnfrac,buffer) 
         call writevar(filen,idies,'burnarea',buffer/1000000.0,istart,icount,ftime,istat)
         end if

! variables which have pft dimension
      
      if (env_debug .gt. 1 ) print *,'writing to '//trim(filenpft)

      istart(3) = 1
      icount(3) = npft
      if ( mlpt .gt. 1 )  then
         istart(4) = 1
         icount(4) = mlpt1
      end if
      istart(ndims+1) = mstep
      icount(ndims+1) = 1
      
      call writevar(filenpft,idiespft,'burnpft',amburnpft,istart,icount,ftime,istat)

      call closefile(idies,istat)
      call closefile(idiespft,istat)

      return
end subroutine wmonthly
