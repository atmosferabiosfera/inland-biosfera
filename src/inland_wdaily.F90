#include "inland_config.h"
#include "inland_compar.h"
#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
subroutine wdaily (jday,nday)
! ---------------------------------------------------------------------
! writes out daily files
!---------------------------------------------------------------
      use inland_parameters
      use inland_control, only: iyear, iyear0, outdir, isimfire, isimveg
      use inland_comatm
      use inland_comfire
      use inland_comwork
      use inland_comsum
      use inland_subgrid
      use inland_comveg
      use inland_comnitr
      use inland_comcrop

      implicit none
!-----------------------------------------------------------------------
! input-output variables
      integer jday,    & ! julian day of the simulation
              nday,    & ! number of days since the start of the simulation
              iyearprev  ! iyear - 1.

! local variables
! parameters used for output
      integer mstep,  & ! this time step for netcdf file
              i,      & ! loop indice 
              j,      & ! loop indice
              idies,  & !
              istat
      integer istart(5), icount(5) ! for writing vars
      real*8 :: tindex(mlpt+1)  ! index used for tiles (last one is average/water)
      real*8 :: pindex(npft)    ! index used for pfts and canopies
      character*8 tmpdate       ! temp. variable to hold date
      character*8 cdate         ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*13 canopies (2) ! canopy definitions
      character*21 tunits       ! units for time
      character*80 dimnames(5)  ! names of dimensions for vars
      character*80 pftdef(npft) ! plant functional type defs (not currently used)
      character*1024 filen        ! file name
      character*4 chyear
      integer ilpt             ! mlpt index
      integer ndims            ! number of dims (including tile, not including pft)
      integer mlpt1            ! number of tiles in output file, including average if mlpt > 1
      real*8 ftime             ! real form of jday


! extra dim variables passed to inifile_dims
      integer, parameter :: numXdim = 2
      integer, parameter :: numXvals = 100 ! pick a big number to make sure we have enough space
      integer nXth(numXdim)
      character*80, dimension(numXDim) :: nameXth, longXth, unitsXth, axisXth
      real*8 valsXth(numXVals,numXdim)

      data istart / 1,1,1,1,1 /, &
           icount / nlon,nlat,1,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub

      iyearprev = iyear - 1

! ---------------------------------------------------------------------
!     current time value, step, time weight
      ftime = nday
      mstep = jday

      chyear = '0000'
      if (iyear .ge. 1000) then
         write(chyear(1:4),'(i4)') iyear
      elseif (iyear .lt. 10) then
         write(chyear(4:4),'(i1)') iyear
      elseif (iyear .lt. 100) then
         write(chyear(3:4),'(i2)') iyear
      else
         write(chyear(2:4),'(i3)') iyear
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

      endif
 
!
! Recording variables in a single file.  
!
      if (myid .eq. 0) then
         filen = trim(outdir)//'/inland-daily-'//chyear//'.nc'
         if (mstep .eq. 1) then
            call inifile(idies,filen,'Inland', &
                         'inland wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         'tile','tile','none','Z',mlpt1,tindex,'none', &
                         tunits,'gregorian',istat)

! inivar
            dimnames(1) = 'longitude'
            dimnames(2) = 'latitude'
            dimnames(3:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(3) = 'tile'

            !gabriel abrahao: Initialize auxlongitude and pids in the output file for reference. 
            if (isparse.eq.1) then
               call inivar(idies,'auxpid','Point ids for sparse mode', &
                       'id',2,dimnames(1:2),istat)
               call inivar(idies,'auxlongitude','Longitudes for sparse mode', &
                       'degrees',2,dimnames(1:2),istat)               
            end if

            call inivar(idies,'rain','average rainfall','mm/day',ndims,dimnames, &
                        istat)
            call inivar(idies,'snow','average snowfall','mm/day',ndims,dimnames, &
                        istat)
!            call inivar(idies,'aet','average aet','mm/day',ndims,dimnames, &
!                        istat)
!            call inivar(idies,'trunoff','average total runoff','mm/day',ndims,   &
!                        dimnames,istat)
            call inivar(idies,'srunoff','average surface runoff','mm/day',ndims, &
                        dimnames,istat)
            call inivar(idies,'drainage','average drainage','mm/day',ndims,      &
                        dimnames,istat)
            call inivar(idies,'tnpptot','average npp','fraction',ndims,dimnames, &
                        istat)
!            call inivar(idies,'wisoi','average soil ice','fraction',ndims,       &
!                        dimnames,istat)            
!            call inivar(idies,'snod','average snow depth','meters',ndims,        &
!                        dimnames,istat)
!            call inivar(idies,'snof','average snow fraction','fraction',ndims,   &
!                         dimnames,istat)
            call inivar(idies,'co2ratio','average co2 ratio','fraction',ndims,   &
                        dimnames,istat)
            call inivar(idies,'co2mic','soil microbe carbon flux','kg/m^2',ndims,&
                        dimnames,istat)
            if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
               call inivar(idies,'Pbio','fire prob due to biomass','fraction',ndims,&
                           dimnames,istat)
               call inivar(idies,'Pmoi','fire prob due to wetness','fraction',ndims,&
                           dimnames,istat)
               call inivar(idies,'Pign','fire prob due to ignition','fraction',ndims,&
                           dimnames,istat)
               call inivar(idies,'Pfire','fire prob','fraction',ndims,&
                           dimnames,istat)
               call inivar(idies,'srate','fire spread rate','km/h',ndims,&
                           dimnames,istat)
               call inivar(idies,'abday','area burned in 1 day','km^2',ndims,&
                           dimnames,istat)
               call inivar(idies,'burnfrac','burned fraction','fraction',ndims,&
                           dimnames,istat)
               call inivar(idies,'burnarea','burned area','km^2',ndims,&
                           dimnames,istat)
               end if
            call inivar(idies,'dnileach','rate of nitrogen leaching', &
                        'kg N ha-1 y-1',ndims,dimnames,istat)

        call closefile(idies,istat)
        
         end if

         call openfile(idies,filen,istat)

      endif

! writevar

      if ( mlpt .gt. 1 )  then
         istart(3) = 1
         icount(3) = mlpt1
      end if
      istart(ndims:) = mstep
      icount(ndims:) = 1

      !gabriel abrahao: Write the pids in the output file for reference.
      if (isparse.eq.1) then
         call writevar(filen,idies,'auxpid',auxpid,istart(1:2),icount(1:2),ftime,istat)
         call writevar(filen,idies,'auxlongitude',auxlonscale,istart(1:2),icount(1:2),ftime,istat)
      end if

      call writevar(filen,idies,'rain',adrain,istart,icount,ftime,istat)
      call writevar(filen,idies,'snow',adsnow,istart,icount,ftime,istat)
! aet
! trunoff
      call writevar(filen,idies,'srunoff',adsrunoff,istart,icount,ftime,istat)
      call writevar(filen,idies,'drainage',addrainage,istart,icount,ftime,istat)
      call writevar(filen,idies,'tnpptot',adtnpptot,istart,icount,ftime,istat)
! wisoi
! snod
! snof
      call writevar(filen,idies,'co2ratio',adco2ratio,istart,icount,ftime,istat)
      call writevar(filen,idies,'co2mic',adco2mic,istart,icount,ftime,istat)
      call writevar(filen,idies,'dnileach',dnileach,istart,icount,ftime,istat)


! fire vars
      if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
         call writevar(filen,idies,'Pbio',adpbio,istart,icount,ftime,istat)
         call writevar(filen,idies,'Pmoi',adpmoi,istart,icount,ftime,istat)
         call writevar(filen,idies,'Pign',adpign,istart,icount,ftime,istat)
         call writevar(filen,idies,'Pfire',adpfire,istart,icount,ftime,istat)
         call writevar(filen,idies,'srate',adsrate,istart,icount,ftime,istat)
         call writevar(filen,idies,'abday',adabday,istart,icount,ftime,istat)
         call writevar(filen,idies,'burnfrac',adburnfrac,istart,icount,ftime,istat)
         call subgrid_calculate_area(adburnfrac,buffer) 
         call writevar(filen,idies,'burnarea',buffer/1000000.0,istart,icount,ftime,istat)
      end if

!    Close file
!    ----------
      call closefile(idies,istat)

!gabriel apagar
write(*,*) 'plai(1,13)',plai(1,13)

!*********************************
! only do if crops are growing
!
   if (isimagro .gt. 0) then

! File created below contains Agro-IBIS variables
      if (myid .eq. 0) then
         filen = 'output/biomass-daily-'//chyear//'.nc'
         if (mstep .eq. 1) then

               ! fill variables passed to inifile_dims
               nameXth(1)  = 'pft'
               longXth(1)  = 'plant functional type'
               unitsXth(1) = 'none'
               axisXth(1)  = 'E'
               nXth(1)     = ecpft-scpft+1!npft
               valsXth(1:ecpft-scpft+1,1) = pindex(scpft:ecpft)
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
              
               call inifile_dims(idies,filen,'daily carbon biomass', &
                         'inland wyearly',cdate,nlonsub,lonscale,nlatsub,latscale, &
                         nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                         tunits,'gregorian',istat)

!            istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
            dimnames(3) = 'pft'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'

            !gabriel abrahao: Initialize auxlongitude and pids in the output file for reference. 
            if (isparse.eq.1) then
               call inivar(idies,'auxpid','Point ids for sparse mode', &
                       'id',2,dimnames(1:2),istat)
               call inivar(idies,'auxlongitude','Longitudes for sparse mode', &
                       'degrees',2,dimnames(1:2),istat)               
            end if

            call inivar(idies,'biomass','biomass status for each pft (kg_C m-2)','kg/m^2', &
                 ndims+1,dimnames,istat)
            call inivar(idies,'cbior','total fine root biomass for each  pft (kg_C m-2)', &
                 'kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'cbiol','total biomass of leaves for each  pft (kg_C m-2)', &
                 'kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'cbios','total stem biomass for each crop pft (kg_C m-2)', &
                 'kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'cbiog','total grain biomass for each crop pft (kg_C m-2)', &
                 'kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'bp','total biomass production for each crop (kg_C m-2)', &
                 'kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'rp','root biomass production for each crop (kg_C m-2)', &
                 '%',ndims+1,dimnames,istat)
            call inivar(idies,'lp','leaf biomass production for each crop (kg_C m-2)', &
                 '%',ndims+1,dimnames,istat)

            call inivar(idies,'lai','leaf area index of each crop (m2_leaf m-2)', &
                 '%',ndims+1,dimnames,istat)            

            call closefile(idies,istat)
         endif

         call openfile(idies,filen,istat)

        istart(3) = 1
        icount(3) = (ecpft-scpft+1) 
        if ( mlpt .gt. 1 )  then
           istart(4) = 1
           icount(4) = mlpt1
        end if
        istart(ndims+1) = mstep
        icount(ndims+1) = 1
 
!        write(*,*) biomass(:,scpft:)

        !gabriel abrahao: Write the pids in the output file for reference.
        if (isparse.eq.1) then
           call writevar(filen,idies,'auxpid',auxpid,istart(1:2),icount(1:2),ftime,istat)
           call writevar(filen,idies,'auxlongitude',auxlonscale,istart(1:2),icount(1:2),ftime,istat)
        end if

        call writevar(filen,idies,'biomass',biomass(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'cbior',cbior(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'cbiol',cbiol(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'cbios',cbios(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'cbiog',cbiog(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'bp',aybprod(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'rp',ayrprod(:,scpft:),istart,icount,ftime,istat)
        call writevar(filen,idies,'lp',aylprod(:,scpft:),istart,icount,ftime,istat)

        call writevar(filen,idies,'lai',plai(:,scpft:),istart,icount,ftime,istat)


      endif

      call closefile(idies,istat)

   endif !check for the crop existence 
!************************** end of AgroIBIS' file biomass.nc **************************


     return
end subroutine wdaily
