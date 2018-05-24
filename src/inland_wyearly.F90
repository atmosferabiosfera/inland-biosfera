#include "inland_config.h"
#include "inland_compar.h"
#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
subroutine wyearly (nday)
! ---------------------------------------------------------------------
! writes out yearly files
! ---------------------------------------------------------------
      use inland_parameters
      use inland_control
      use inland_comfire
      use inland_comcrop
      use inland_comveg
      use inland_comwork
      use inland_comnitr
      use inland_comsum
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
      integer :: idies, idiespft, idiestile 
      integer :: mstep, istat, i, k, zero
      integer :: istart(5), icount(5)  ! for writing vars
      real*8 :: pindex(npft)    ! index used for pfts and canopies
      real*8 :: tindex(mlpt+1)  ! index used for tiles (last one is average/water)
      character*8 tmpdate       ! temp. variable to hold date
      character*8 cdate         ! date to use in history attribute in files
      character*13 canopies(2)  ! canopy definitions
      character*1024 canopies2  ! canopy definitions
      character*21 tunits       ! time units
      character*80 dimnames(5)  ! names of dimensions for vars
      character*60 pftdef(npft) ! plant functional type defs
      character*1024 pftdef2
      character*60 vegtypedef(15) ! vegtype defs (15 currently)
      character*1024 vegtypedef2
      character*1024 tiledef
      character*1024 filen, filenpft, filentile ! file name
      character*4 chyear
      integer ilpt             ! mlpt index
      integer ndims          ! number of dims (including tile, not including pft)
      integer mlpt1            ! number of tiles in output file, including average if mlpt > 1
      real*8 ftime         ! real form of nday
      character*1 endline   ! endline character

! extra dim variables passed to inifile_dims
      integer, parameter :: numXdim = 2
      integer, parameter :: numXvals = 100 ! pick a big number to make sure we have enough space
      integer nXth(numXdim)
      character*80, dimension(numXDim) :: nameXth, longXth, unitsXth, axisXth
      real*8 valsXth(numXVals,numXdim)

      data istart / 1,1,1,1,1 /,  &
           icount / nlon,nlat,1,1,1 /
! ---------------------------------------------------------------------
! FIXME: no point in using 'data' above if we change it right away to a
!       variable that is not compatible with 'data' assignment.
      icount(1) = nlonsub
      icount(2) = nlatsub

      zero = 0 ! you kiddin!
      endline = char(10) ! this is not read properly by ncdump<4.3 when using netcdf-4 files

! current time value and step: make ftime Jan 1 of this year
! instead of Dec 31
      ftime = nday - ndaypy + 1
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

! how many dims?
      if ( mlpt .gt. 1 ) then
         ndims = 4
         mlpt1 = mlpt + 1
      else
         ndims = 3
         mlpt1 = 1
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

! define plant functional types
         do i = 1, npft
            pindex(i) = i
         end do

! define tile dimensions
         do ilpt = 1, mlpt1
            tindex(ilpt) = ilpt
         end do
         ilpt = 1

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

      end if ! (mstep .eq. 1)

      filen = trim(outdir)//'/inland-yearly-'//chyear//'.nc'
      filenpft = trim(outdir)//'/inland-yearly-pft-'//chyear//'.nc'

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
                         'inland wyearly',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                         tunits,'gregorian',istat)
               call inifile_dims(idiespft,filenpft,'Inland', &
                         'inland wyearly pft',cdate,nlonsub,lonscale,nlatsub,latscale, &
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

            call inivar(idies,'npptot','total npp','kg m-2 year-1',ndims,dimnames, &
                        istat)
            call inivar(idies,'anpptot','total above-ground npp', &
                        'kg m-2 year-1',ndims,dimnames,istat)
            call inivar(idies,'aet','average evapotranspiration','mm/year',ndims, &
                        dimnames,istat)
            call inivar(idies,'trunoff','total runoff','mm/year',ndims,dimnames,  &
                        istat)
            call inivar(idies,'srunoff','surface runoff','mm/year',ndims,dimnames, &
                        istat)
            call inivar(idies,'drainage','drainage','mm/year',ndims,dimnames,&
                        istat)
            call inivar(idies,'wsoi','average soil moisture','fraction',ndims,     &
                        dimnames,istat)
            call inivar(idies,'wisoi','average soil ice','fraction',ndims,dimnames,&
                        istat)
            call inivar(idies,'vwc','average volumetric water content', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'awc','average volumetric water content','cm',ndims, &
                        dimnames,istat)
            call inivar(idies,'tsoi','average soil temperature','degrees C',ndims, &
                        dimnames,istat)
            call inivar(idies,'totlaiu','total lai for upper canopy', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'totlail','total lai for lower canopy', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'totbiou','total biomass for upper canopy', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totbiol','total biomass for lower canopy', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totfall','total litterfall','kg/m^2',ndims,dimnames,&
                        istat)
            call inivar(idies,'rootbio','total live root biomass carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totalit','total above ground litter carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totrlit','total below ground litter carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totcsoi','total soil carbon w/o litter', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totcmic','total microbial carbon','kg/m^2',ndims, &
                        dimnames,istat)
            call inivar(idies,'totanlit','total above ground litter nitrogen', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totrnlit','total below ground litter nitrogen', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'totnsoi','total soil nitrogen w/o litter', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'nmintot','total nitrogen mineralization', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'neetot','total net ecosystem echange carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'caccount','end of year carbon correction', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'co2mic','total microbe respiration carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'co2root','total root respiration carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'co2soi','total soil respiration carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'fu','fractional cover of upper canopy', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'fl','fractional cover of lower canopy', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'gdd8clim','climate gdd8', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd10clim','climate gdd10', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd12clim','climate gdd12', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd0clim','climate gdd0', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'currgdd8','current year gdd8', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'currgdd10','current year gdd10', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'currgdd12','current year gdd12', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'currgdd0','current year gdd0', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gddfzcorn','current year gdd between freeze events', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gddfzsoy','current year gdd between freeze events', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gddfzsgc','current year gdd between freeze events', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gsdays','number of days between freeze events', &
                        'days',ndims,dimnames,istat)
            call inivar(idies,'iniday','last day in spring for freeze', &
                        'day',ndims,dimnames,istat)
            call inivar(idies,'falll', 'total annual leaf litterfall', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'fallr', 'total below ground annual root turnover', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'fallw','total wood litterfall', &
                       'kg/m^2',ndims,dimnames,istat)

            call inivar(idies,'nimmob','total nitrogen immobilized', &
                       'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'nfixnat','nitrogen fixation by natural vegetation', &
                       'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'ftot','annual total inorganic nitrogen leached from soil profile (kg n  m-2 y-1)', &
                       'kg/hectare',ndims,dimnames,istat)
            call inivar(idies,'no3leach','annual total nitrate leaching', &
                       'kg/hectare',ndims,dimnames,istat)
            call inivar(idies,'nbalance','annual inorganic nitrogen balance', &
                       'kg/hectare',ndims,dimnames,istat)
            call inivar(idies,'natvegnup','total annual inorganic nitrogen uptake by natural veg', &
                       'kg/m2/yr',ndims,dimnames,istat)


           call inivar(idies,'csoislo','total soil carbon in slow pool', &
                      'kg/m^2',ndims,dimnames,istat)
           call inivar(idies,'csoipas','total soil carbon in passive pool', &
                      'kg/m^2',ndims,dimnames,istat)

           call inivar_byte(idies,'vegtype0','vegetation type','none',&
                            ndims,dimnames,istat)

!            call inivar(idies,'vegtype0','vegetation type','none',ndims,dimnames, &
!                        istat)
            call inivar_byte(idies,'vegtype','vegetation type - static','none',&
                            ndims,dimnames,istat)
            call inivar_byte(idies,'landusetype','land use (1-4)','none',ndims,dimnames, &
                            istat)
            if ( mlpt .gt. 1 ) then
               call inivar(idies,'tilefrac','tile fraction','none',ndims,dimnames, &
                    istat)
            end if
            if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
               call inivar(idies,'Pfire','Fire prob', &
                           'fraction',ndims,dimnames,istat)
               call inivar(idies,'srate','fire spread rate', &
                           'km/h',ndims,dimnames,istat)
               call inivar(idies,'vegfrac','vegetated fraction', &
                           'km^2',ndims,dimnames,istat)
               call inivar(idies,'burnfrac','burned fraction', &
                           'km^2',ndims,dimnames,istat)
               call inivar(idies,'burnarea','burned area', &
                           'km^2',ndims,dimnames,istat)
            end if
            call inivar_byte(idies,'croptype','crop type', 'none',ndims,dimnames,istat)

! variables which have pft dimension
! for reading these files with cdo, when tiling is used, see "Subgrid Tiling" in docs/README.advanced

            dimnames(3) = 'pft'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'
            call inivar_byte(idiespft,'exist','existence for each pft','none',ndims+1, &
                        dimnames,istat)
!            call inivar(idiespft,'exist','existence for each pft','none',ndims+1, &
!                        dimnames,istat)
            call inivar(idiespft,'plai','leaf area index for each pft','fraction',&
                        ndims+1,dimnames,istat)
            call inivar(idiespft,'biomass','biomass for each pft','kg/m^2',ndims+1, &
                        dimnames,istat)
            call inivar(idiespft,'npp','npp of carbon for each pft', &
                        'kg m-2 year-1',ndims+1,dimnames,istat)
            call inivar(idiespft,'burnpft','burned fraction for each pft', &
                        'km^2',ndims+1,dimnames,istat)

            call closefile(idies,istat)
            call closefile(idiespft,istat)

         end if

         call openfile(idies,filen,istat)
         call openfile(idiespft,filenpft,istat)


      endif

      if (env_debug .gt. 1 ) print *,'writing to '//trim(filen)

! these vars have no subgrid average/total
! vegetation type
      istart(3) = 1
      icount(3) = mlpt
      istart(4) = mstep
      icount(4) = 1
      call writevar(filen,idies,'vegtype0',vegtype0,istart,icount,ftime,istat)
      call writevar(filen,idies,'vegtype',xinveg,istart,icount,ftime,istat)

! land use
      call writevar(filen,idies,'landusetype',landusetype,istart,icount,ftime,istat)
      if ( mlpt .gt. 1 ) call writevar(filen,idies,'tilefrac',tilefrac,istart,icount,ftime,istat)

! last tilefrac tile is waterfrac
      if ( mlpt .gt. 1 ) then
         istart(3) = mlpt+1
         icount(3) = 1
         istart(4) = mstep
         icount(4) = 1
         
         call writevar(filen,idies,'tilefrac',waterfrac,istart,icount,ftime,istat)
      endif

! these vars have subgrid average
      if ( mlpt .gt. 1 )  then
         istart(3) = 1
         icount(3) = mlpt1
      end if
      istart(ndims) = mstep
      icount(ndims) = 1

      call writevar(filen,idies,'npptot',aynpptot,istart,icount,ftime,istat)
      call writevar(filen,idies,'anpptot',ayanpptot,istart,icount,ftime,istat)

! evapotranspiration

      call writevar(filen,idies,'aet',ayaet,istart,icount,ftime,istat)

! trunoff, srunoff, drainage, rratio, tratio
      call writevar(filen,idies,'trunoff',aytrunoff,istart,icount,ftime,istat)
      call writevar(filen,idies,'srunoff',aysrunoff,istart,icount,ftime,istat)
      call writevar(filen,idies,'drainage',aydrainage,istart,icount,ftime,istat)

! soil moisture, soil ice, volumetric water content, plant available water
      call writevar(filen,idies,'wsoi',aywsoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'wisoi',aywisoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'vwc',ayvwc,istart,icount,ftime,istat)
      call writevar(filen,idies,'awc',ayawc,istart,icount,ftime,istat)

! soil temperature
      call writevar(filen,idies,'tsoi',aytsoi,istart,icount,ftime,istat)

! solar radiation
! albedo
! upward and downward infrared radiation
! sensible heat flux
! latent heat flux

! lai, by pft, total upper canopy, total lower canopy
      call writevar(filen,idies,'totlaiu',totlaiu,istart,icount,ftime,istat)
      call writevar(filen,idies,'totlail',totlail,istart,icount,ftime,istat)

! biomass, by pft, upper canopy, lower canopy
      call writevar(filen,idies,'totbiou',totbiou,istart,icount,ftime,istat)
      call writevar(filen,idies,'totbiol',totbiol,istart,icount,ftime,istat)

! soil carbon: rootbio, totalit, totrlit, totcsoi, totcmic
      call writevar(filen,idies,'rootbio',ayrootbio,istart,icount,ftime,istat)
      call writevar(filen,idies,'totalit',ayalit,istart,icount,ftime,istat)
      call writevar(filen,idies,'totrlit',ayblit,istart,icount,ftime,istat)
      call writevar(filen,idies,'totcsoi',aycsoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'totcmic',aycmic,istart,icount,ftime,istat)

! soil nitrogen: totanlit, totrnlit, totnsoi, nminott
      call writevar(filen,idies,'totanlit',ayanlit,istart,icount,ftime,istat)
      call writevar(filen,idies,'totrnlit',aybnlit,istart,icount,ftime,istat)
      call writevar(filen,idies,'totnsoi',aynsoi,istart,icount,ftime,istat)
      call writevar(filen,idies,'nmintot',aynmintot,istart,icount,ftime,istat)

! total litter
! total litterfall
      call writevar(filen,idies,'totfall',totfall,istart,icount,ftime,istat)

! total soil carbon in slow pool
! total soil carbon in passive pool
! co2 carbon exchange: net ecosystem, microbial resp, root resp, soil resp
      call writevar(filen,idies,'neetot',ayneetot,istart,icount,ftime,istat)
      call writevar(filen,idies,'caccount',caccount,istart,icount,ftime,istat)
      call writevar(filen,idies,'co2mic',ayco2mic,istart,icount,ftime,istat)
      call writevar(filen,idies,'co2root',ayco2root,istart,icount,ftime,istat)
      call writevar(filen,idies,'co2soi',ayco2soi,istart,icount,ftime,istat)

! fire disturbance regime

! fractional cover of upper and lower canopies
      call writevar(filen,idies,'fu',fu,istart,icount,ftime,istat)
      call writevar(filen,idies,'fl',fl,istart,icount,ftime,istat)
      call writevar(filen,idies,'gdd8clim',gdd8,istart,icount,ftime,istat)
      call writevar(filen,idies,'gdd10clim',gdd10,istart,icount,ftime,istat)
      call writevar(filen,idies,'gdd12clim',gdd12,istart,icount,ftime,istat)
      call writevar(filen,idies,'gdd0clim',gdd0c,istart,icount,ftime,istat)
      call writevar(filen,idies,'currgdd8',gdd8this,istart,icount,ftime,istat)
      call writevar(filen,idies,'currgdd10',gdd10this,istart,icount,ftime,istat)
      call writevar(filen,idies,'currgdd12',gdd12this,istart,icount,ftime,istat)
      call writevar(filen,idies,'currgdd0',gdd0cthis,istart,icount,ftime,istat)

      call writevar(filen,idies,'gddfzcorn',gddcorn,istart,icount,ftime,istat)
      call writevar(filen,idies,'gddfzsoy',gddsoy,istart,icount,ftime,istat)
      call writevar(filen,idies,'gddfzsgc',gddsgc,istart,icount,ftime,istat)
      call writevar(filen,idies,'gsdays',gsdays,istart,icount,ftime,istat)
      buffer = iniday
      call writevar(filen,idies,'iniday',buffer,istart,icount,ftime,istat)

      call writevar(filen,idies,'falll',falll,istart,icount,ftime,istat)
      call writevar(filen,idies,'fallr',fallr,istart,icount,ftime,istat)
      call writevar(filen,idies,'fallw',fallw,istart,icount,ftime,istat)

      call writevar(filen,idies,'nimmob',ayimmtot,istart,icount,ftime,istat)
      call writevar(filen,idies,'nfixnat',yfixsoin,istart,icount,ftime,istat)
      call writevar(filen,idies,'ftot',ftot,istart,icount,ftime,istat)
      call writevar(filen,idies,'no3leach',yno3leach,istart,icount,ftime,istat)
      call writevar(filen,idies,'nbalance',snbalance,istart,icount,ftime,istat)
      call writevar(filen,idies,'natvegnup',totnvegn,istart,icount,ftime,istat)
      call writevar(filen,idies,'csoislo',csoislo,istart,icount,ftime,istat)
      call writevar(filen,idies,'csoipas',csoipas,istart,icount,ftime,istat)

! vegetation type
      istart(3) = 1
      icount(3) = mlpt
      istart(4) = mstep
      icount(4) = 1
      call writevar(filen,idies,'vegtype0',vegtype0,istart,icount,ftime,istat)

      call writevar(filen,idies,'vegtype',xinveg,istart,icount,ftime,istat)

! land use
      call writevar(filen,idies,'landusetype',landusetype,istart,icount,ftime,istat)
      call writevar(filen,idies,'croptype',croptype,istart,icount,ftime,istat)


! sapwood fraction
! bottom and top of vegetation canopies

! fire disturbance regime
      if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
         call writevar(filen,idies,'Pfire',aypfire,istart,icount,ftime,istat)
         call writevar(filen,idies,'srate',aysrate,istart,icount,ftime,istat)
         call writevar(filen,idies,'vegfrac',vegfrac,istart,icount,ftime,istat)
         call writevar(filen,idies,'burnfrac',ayburnfrac,istart,icount,ftime,istat)
         !call writevar(filen,idies,'burnarea',ayburnfrac*garea/1000000.0,&
         !              istart,icount,ftime,istat)
         call subgrid_calculate_area(ayburnfrac,buffer) 
         call writevar(filen,idies,'burnarea',buffer/1000000.0,istart,icount,ftime,istat)
      end if


! variables which have pft dimension
      
      if (env_debug .gt. 1 ) print *,'writing to '//trim(filenpft)

      istart(3) = 1
      icount(3) = npft
      if ( mlpt .gt. 1 )  then
         istart(4) = 1
         icount(4) = mlpt
      end if
      istart(ndims+1) = mstep
      icount(ndims+1) = 1
      
      call writevar(filenpft,idiespft,'exist',exist,istart,icount,ftime,istat)
      icount(4) = mlpt1
      call writevar(filenpft,idiespft,'plai',plai,istart,icount,ftime,istat)
      call writevar(filenpft,idiespft,'biomass',biomass,istart,icount,ftime,istat)
      call writevar(filenpft,idiespft,'npp',aynpp,istart,icount,ftime,istat)
      call writevar(filenpft,idiespft,'burnpft',ayburnpft,istart,icount,ftime,istat)

      call closefile(idies,istat)
      call closefile(idiespft,istat)


! -----------------------------------------------------------------------
! 1 file with all tile vars, by tile level
      
      if ( ( mlpt .gt. 1 ) .or. ( trim(hrmapfile) .ne. '' ) ) then

      filentile = trim(outdir)//'/inland-tiles-'//chyear//'.nc'
      !filen4 = trim(outdir)//'/inland-tiles-nofill-'//chyear//'.nc'

      if (env_debug .gt. 1 ) print *,'creating '//trim(filentile)

      if (myid .eq. 0) then
         if (mstep .eq. 1) then

            ! create file to store tile/veg information, so main file is still compat.
            ! eventually merge into main file
            call inifile(idiestile,filentile,'Inland', &
                         'inland wyearly - tiles',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         'tile','tile','none','Z',mlpt+1,tindex,'none', &
                         tunits,'gregorian',istat)

            dimnames(3) = 'tile'
            dimnames(4) = 'time'
            call inivar_byte(idiestile,'vegtype0','vegetation type (dynamic)','none',&
                        ndims,dimnames,istat)
            call inivar_byte(idiestile,'vegtype','vegetation type (static)','none',&
                        ndims,dimnames,istat)
            call inivar(idiestile,'tilefrac','tile fraction','none',ndims,dimnames, &
                         istat)
            call inivar_byte(idiestile,'landusetype','land use (1-4)','none',ndims,dimnames, &
                            istat)
            call inivar_int(idiestile,'itilechild','tile child index','none', &
                           ndims,dimnames,istat)
            call inivar_int(idiestile,'itileparent','tile parent index','none', &
                           ndims,dimnames,istat)
            dimnames(3) = 'time'
            dimnames(4) = 'none'
            call inivar_byte(idiestile,'ntilechild','number of children tiles','none',3,dimnames, &
                            istat)
            ! ! this file is just for debugging purposes when reading hrmap
            ! if ( trim(hrmapfile) .ne. '' ) then
            ! call inivarbyte(idiesndims,'vegtype','vegetation type','none',&
            !                 ndims,dimnames,istat)
            ! call inivar(idiesndims,'tilefrac','tile fraction','none',ndims,dimnames, &
            !              istat)
            ! end if
            !dimnames(3) = 'time'

            call closefile(idiestile,istat)

         end if

         call openfile(idiestile,filentile,istat)

      end if

      istart(3) = 1
      icount(3) = mlpt
      istart(4) = mstep
      icount(4) = 1

! TODO fix when input data is integer, not double!!!
! vegetation type - dynamic
      call writevar(filentile,idiestile,'vegtype0',vegtype0,istart,icount,ftime,istat)

! vegetation type - static
      call writevar(filentile,idiestile,'vegtype',xinveg,istart,icount,ftime,istat)

! subgrid tile fraction
      call writevar(filentile,idiestile,'tilefrac',tilefrac,istart,icount,ftime,istat)

! land use
      call writevar(filen,idiestile,'landusetype',landusetype,istart,icount,ftime,istat)


! itileparent
      buffer = itileparent ! need to copy integer to real variable
      call writevar(filentile,idiestile,'itileparent',buffer,istart,icount,ftime,istat)

! itilechild
      ! itilechild is a 2D var, need to loop over mlpt
      do ilpt = 1, mlpt

         istart(3) = ilpt
         icount(3) = 1
         istart(4) = mstep
         icount(4) = 1

         buffer1 = itilechild(:,ilpt)
         call writevar(filentile,idiestile,'itilechild',buffer1,istart,icount,ftime,istat)

      end do

! ntilechild - 1:npoi1
      istart(3) = mstep
      icount(3) = 1
      istart(4) = 1
      icount(4) = 1

      buffer1 = ntilechild
      call writevar(filentile,idiestile,'ntilechild',buffer1,istart,icount,ftime,istat)


! last tile is waterfrac (vegtype=30)
      istart(3) = mlpt+1
      icount(3) = 1
      istart(4) = mstep
      icount(4) = 1

      buffer1 = 30.
      !print *,cdummy
      call writevar(filentile,idiestile,'vegtype',buffer1,istart,icount,ftime,istat)
      call writevar(filentile,idiestile,'tilefrac',waterfrac,istart,icount,ftime,istat)


      ! close file
      call closefile(idiestile,istat)

      end if !( ( mlpt .gt. 1 ) .or. ( trim(hrmapfile) .ne. '' ) )

!done
!


!*********************crop yearly**********************************
!*********************************
! only do if crops are growing
!
   if (isimagro .gt. 0) then


     if (myid .eq. 0) then
         filen = 'output/crop-yearly-'//chyear//'.nc'
         if (mstep .eq. 1) then

               ! fill variables passed to inifile_dims
               nXth(1)     = ecpft-scpft+1!npft
               valsXth(1:ecpft-scpft+1,1) = pindex(scpft:ecpft)
               ! tile already filled

               call inifile_dims(idies,filen,'annual crop growth variables', &
                         'inland wyearly',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                         tunits,'gregorian',istat)

            istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
            
            dimnames(3) = 'pft'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'

         call inivar(idies,'cropy','crop cycle number in a cycle (1s planting - last ratoon)', &
                    'day of year',ndims+1,dimnames,istat)

         call inivar(idies,'plantdate','crop planting date -day', &
                    'day of year',ndims+1,dimnames,istat)

         call inivar(idies,'hdate','crop harvest date - day', &
                    'day of year',ndims+1,dimnames,istat)

         call inivar(idies,'dpp',' number of days of crop cycle - day', &
                    'number of days',ndims+1,dimnames,istat)

         call inivar(idies,'harvidx','harvest index - fraction', &
                    'fraction',ndims+1,dimnames,istat)

         call inivar(idies,'croplaimx','maximum crop lai- m^2/m^2', &
                    'm^2/m^2',ndims+1,dimnames,istat)

         call inivar(idies,'grainn','nitrogen removed by grain - kg/ha', &
                    'kg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'cropyld','crop yield - grain(t/ha) or stalk (t/ha)', &
                    't/ha',ndims+1,dimnames,istat)

         call inivar(idies,'dmyield','crop yield dry matter - Mg/ha', &
                    'Mg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'dmleaf','leaf dry matter - Mg/ha', &
                    'Mg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'dmstem','stem dry matter - Mg/ha', &
                    'Mg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'dmroot','root dry matter - Mg/ha', &
                    'Mg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'dmresidue','aboveground residue dry matter -Mg/ha', &
                    'Mg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'dmcrop','total crop dry matter - Mg/ha', &
                    'Mg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'residuen','total nitrogen in residue - kg/ha', &
                    'kg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'nconcl','leaf nitrogen concentration - %', &
                    'percent',ndims+1,dimnames,istat)

         call inivar(idies,'nconcs','stem nitrogen concentration - %', &
                    'percent',ndims+1,dimnames,istat)

         call inivar(idies,'nconcr','root nitrogen concentration - %', &
                    'percent',ndims+1,dimnames,istat)

         call inivar(idies,'nconcg','grain nitrogen concentration - %', &
                    'percent',ndims+1,dimnames,istat)

         call inivar(idies,'cropn','total crop nitrogen uptake - kg/ha', &
                    'kg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'cropfixn','nitrogen fixation - kg/ha', &
                    'kg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'cntops','cn ratio of aboveground residue ', &
                    'dimensionless',ndims+1,dimnames,istat)

         call inivar(idies,'cnroot','cn ratio of  fineroots ', &
                    'dimensionless',ndims+1,dimnames,istat)

         call inivar(idies,'fertilizer','fertilizer applied - kg/ha', &
                    'kg/ha',ndims+1,dimnames,istat)

         call inivar(idies,'grainday','grainday - day', &
                    'day',ndims+1,dimnames,istat)

!         call inivar(idies,'hybrid','hybridgdd - degrees C', &
!                    'degrees C',ndims+1,dimnames,istat)

         call inivar(idies,'gddplant','accumulated growing degree days between planting and harvest', &
                    'degrees C',ndims+1,dimnames,istat)
 
         call inivar(idies,'crmclim','CRM rating for climatology -days', &
                    'days',ndims+1,dimnames,istat)

         call inivar(idies,'crmact','best CRM rating for this year -days', &
                    'days',ndims+1,dimnames,istat)

         call inivar(idies,'crmplant','CRM rating for this year based on planting date -d', &
                    'days',ndims+1,dimnames,istat)

         call closefile(idies,istat)
       endif
         call openfile(idies,filen,istat)
      endif

!      istart(3) = scpft
!      icount(3) = (ecpft-scpft+1)
      istart(3) = 1
      icount(3) = (ecpft-scpft+1)
      if ( mlpt .gt. 1 )  then
         istart(4) = 1
         icount(4) = mlpt1
      end if
      istart(ndims+1) = mstep
      icount(ndims+1) = 1

      call writevar(filen,idies,'cropy',cropout(:,scpft:,29),istart,icount,ftime,istat)
      call writevar(filen,idies,'plantdate',cropout(:,scpft:,1),istart,icount,ftime,istat)
      call writevar(filen,idies,'hdate',cropout(:,scpft:,2),istart,icount,ftime,istat)
      call writevar(filen,idies,'dpp',cropout(:,scpft:,28),istart,icount,ftime,istat)
      call writevar(filen,idies,'harvidx',cropout(:,scpft:,3),istart,icount,ftime,istat)
      call writevar(filen,idies,'croplaimx',cropout(:,scpft:,4),istart,icount,ftime,istat)
      call writevar(filen,idies,'grainn',cropout(:,scpft:,5),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropyld',cropout(:,scpft:,6),istart,icount,ftime,istat)
      call writevar(filen,idies,'dmyield',cropout(:,scpft:,7),istart,icount,ftime,istat)
      call writevar(filen,idies,'dmleaf',cropout(:,scpft:,8),istart,icount,ftime,istat)
      call writevar(filen,idies,'dmstem',cropout(:,scpft:,9),istart,icount,ftime,istat)
      call writevar(filen,idies,'dmroot',cropout(:,scpft:,10),istart,icount,ftime,istat)
      call writevar(filen,idies,'dmresidue',cropout(:,scpft:,11),istart,icount,ftime,istat)
      call writevar(filen,idies,'dmcrop',cropout(:,scpft:,12),istart,icount,ftime,istat)
      call writevar(filen,idies,'residuen',cropout(:,scpft:,13),istart,icount,ftime,istat)
      call writevar(filen,idies,'nconcl',cropout(:,scpft:,14),istart,icount,ftime,istat)
      call writevar(filen,idies,'nconcs',cropout(:,scpft:,15),istart,icount,ftime,istat)
      call writevar(filen,idies,'nconcr',cropout(:,scpft:,16),istart,icount,ftime,istat)
      call writevar(filen,idies,'nconcg',cropout(:,scpft:,13),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropn',cropout(:,scpft:,18),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropfixn',cropout(:,scpft:,19),istart,icount,ftime,istat)
      call writevar(filen,idies,'cntops',cropout(:,scpft:,20),istart,icount,ftime,istat)
      call writevar(filen,idies,'cnroot',cropout(:,scpft:,21),istart,icount,ftime,istat)
      call writevar(filen,idies,'fertilizer',cropout(:,scpft:,22),istart,icount,ftime,istat)
      call writevar(filen,idies,'grainday',cropout(:,scpft:,27),istart,icount,ftime,istat)
      call writevar(filen,idies,'gddplant',cropout(:,scpft:,30),istart,icount,ftime,istat)
      call writevar(filen,idies,'crmclim',cropout(:,scpft:,24),istart,icount,ftime,istat)
      call writevar(filen,idies,'crmact',cropout(:,scpft:,25),istart,icount,ftime,istat)
      call writevar(filen,idies,'crmplant',cropout(:,scpft:,26),istart,icount,ftime,istat)

      call closefile(idies,istat)

endif !check for the crop existence


#if 0
!**********zbot and ztop *****************
!*********************************
! only do if crops are not growing
!
   if (isimagro .gt. 0) then

               nameXth(1)  = 'canopy'
               longXth(1)  = 'canopy'
               unitsXth(1) = 'none'
               axisXth(1)  = 'E'
               nXth(1)     = 2
               valsXth(1:npft,1) = pindex

     if (myid .eq. 0) then
      filen = 'output/zcanopy'//chyear//'.nc'
      if (mstep .eq. 1) then
         call inifile_dims(idies,filen,'annual height of vegetation canopies', &
                         'inland wyearly',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                         tunits,'gregorian',istat)
! add global attribute to define canopies with text, use netcdf low-level com
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def', &
         26+5,'1='//canopies(1)//' 2='//canopies(2))

         dimnames(3) = 'canopy'
         dimnames(4:) = 'time'
         if ( mlpt .gt. 1 ) dimnames(4) = 'tile'

           call inivar(idies,'zbot','bottom heights of lower and upper canopies', &
                      'meters',ndims+1,dimnames,istat)
           call inivar(idies,'ztop','top heights of lower and upper canopies', &
                      'meters',ndims+1,dimnames,istat)
            
         call closefile(idies,istat)
       endif
         call openfile(idies,filen,istat)

      istart(3) = 1
      icount(3) = 2
      if ( mlpt .gt. 1 )  then
         istart(4) = 1
         icount(4) = mlpt1
      end if
      istart(ndims+1) = mstep
      icount(ndims+1) = 1


   call writevar(filen,idies,'zbot',buffer,istart,icount,ftime,istat)
   call writevar(filen,idies,'ztop',ztop(1:,:),istart,icount,ftime,istat)

      if ( mlpt .gt. 1 )  then
         istart(4) = 1
         icount(4) = mlpt1
      end if
      istart(ndims+1) = mstep
      icount(ndims+1) = 1



do i=1,2
      istart(3) = i
      icount(3) = 1
buffer(:) = zbot(1,1)

   call writevar(filen,idies,'zbot',buffer,istart,icount,ftime,istat)  
   call writevar(filen,idies,'ztop',ztop,istart,icount,ftime,istat)

end do

         call closefile(idies,istat)
  endif

endif !check for the crop existence
#endif
      return
end subroutine wyearly
