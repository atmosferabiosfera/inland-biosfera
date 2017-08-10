#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine gdiag(snorth, ssouth, swest, seast, nrun, iyrlast, co2init, o2init, soilcspin)
      use inland_parameters
      use inland_control, only: iyear, iyear0, isimveg, isimfire, isimco2, &
                                indir, infile, datadir, outdir, isimland, &
                                isinfilt, isimrwu, itauw,isla,ica,ivmax
      use inland_comveg
      use inland_comsum
      use inland_subgrid
      use inland_comnitr
      use inland_comcrop, only: isimagro, icroptype

      implicit none

! ------------------------------Variables--------------------------------
! local variables
      integer i, j, ilpt

      real*8 :: gnee,      &  ! domain total nee (gt-c/yr)
                gnpp,      &  ! domain total npp (gt-c/yr)
                ggpp,      &  !  domain total gpp (gt-c/yr)
                gbiomass,  &  ! domain total biomass (gt-c)
                galitc,    &  !  domain total aboveground litter carbon (gt-c)
                gblitc,    &  ! domain total belowground litter carbon (gt-c)
                gsoic,     &  ! domain total soil carbon (gt-c)
                gco2soi,   &  ! domain total soil surface co2 flux (gt-c)
                galitn,    &  ! domain total aboveground litter nitrogen (gt-c)
                gblitn,    &  ! domain total belowground litter nitrogen (gt-c)
                gsoin,     &  ! domain total soil nitrogen (gt-c)
                gprcp,     &  ! domain average annual precipitation (mm/yr)
                gaet,      &  ! domain average annual evapotranspiration (mm/yr)
                gt,        &  ! domain average annual transpiration (mm/yr)
                gtrunoff,  &  ! domain average total runoff (mm/yr)
                gsrunoff,  &  ! domain average surface runoff (mm/yr)
                gdrainage, &  ! domain average drainage (mm/yr)
                gdwtot,    &  !   "      "     water recharge (mm/yr)
                gco2mic,   &  !
                gcdist,    &  !
                gtotfall,  &  !
                gtarea,    &  ! total land area of the domain (m**2)
                aratio,    &  ! aet / prcp ratio
                rratio,    &  ! runoff / prcp ratio
                sratio,    &  ! surface runoff / drainage ratio
                tratio,    &  ! transpiration / aet ratio
                gleach        ! domain total leaching rate
      
      real*8 scaling          ! scaling to apply to each point = garea(i) * tilefrac(i)

      real*8  snorth,      & ! north latitude for subsetting std grid
              ssouth,      & ! south latitude for subsetting std grid
              swest,       & ! west longitude for subsetting std grid
              seast,       & ! east longitude for subsetting std grid
              co2init,     & ! initial co2 concentration in mol/mol
              o2init         ! initial o2 concentration in mol/mol

      integer iyrlast,     & ! last year of previous run (for restart)
              nrun,        & ! # of years in this run
              soilcspin      ! 0: no spinup procedure for soil c  1: acceleration procedure used

! initialize variables
      gtarea    = 0.0
      gnee      = 0.0
      gnpp      = 0.0
      ggpp      = 0.0
      gbiomass  = 0.0
      galitc    = 0.0
      gblitc    = 0.0
      gsoic     = 0.0
      gco2soi   = 0.0
      galitn    = 0.0
      gblitn    = 0.0
      gsoin     = 0.0
      gprcp     = 0.0
      gaet      = 0.0
      gt        = 0.0
      gtrunoff  = 0.0
      gsrunoff  = 0.0
      gdrainage = 0.0
      gdwtot    = 0.0
      gleach    = 0.0
      gco2mic   = 0.0
      gcdist    = 0.0
      gtotfall  = 0.0

! calculate domain-wide totals and averages
      do j = 1, npoi1
         do ilpt = 1, mlpt

            i = subgrid_get_index(j,ilpt)
            ! for now garea is the area of the entire cell, if this changes must remove tilefrac here
            scaling = garea(i) * tilefrac(i)

            ayrratio(i) = min (dble(1.0), max (dble(0.0), aytrunoff(i)) / max (dble(0.1), ayprcp(i)))
            aytratio(i) = min (dble(1.0), max (dble(0.0), aytrans(i))  / max (dble(0.1), ayaet(i)))

            gtarea    = gtarea    + scaling 
            gnee      = gnee      + scaling * ayneetot(i) * 1.e-12
            gnpp      = gnpp      + scaling * aynpptot(i) * 1.e-12
            ggpp      = ggpp      + scaling * aygpptot(i) * 1.e-12
            gbiomass  = gbiomass  + scaling * (totbiou(i) + totbiol(i)) * 1.e-12
            galitc    = galitc    + scaling * totalit(i)  * 1.e-12
            gblitc    = gblitc    + scaling * totrlit(i)  * 1.e-12
            gsoic     = gsoic     + scaling * totcsoi(i)  * 1.e-12
            gco2soi   = gco2soi   + scaling * ayco2soi(i) * 1.e-12
            galitn    = galitn    + scaling * ayanlit(i)  * 1.e-12
            gblitn    = gblitn    + scaling * aybnlit(i)  * 1.e-12
            gsoin     = gsoin     + scaling * aynsoi(i)   * 1.e-12
            gprcp     = gprcp     + scaling * ayprcp(i)
            gaet      = gaet      + scaling * ayaet(i)
            gt        = gt        + scaling * aytrans(i)
            gsrunoff  = gsrunoff  + scaling * aysrunoff(i)
            gdrainage = gdrainage + scaling * aydrainage(i)
            gdwtot    = gdwtot    + scaling * dwtot(i)    * nytime(1)
            gleach    = gleach    + scaling * yrleach(i)  * 1.e-12
            gco2mic   = gco2mic   + scaling * ayco2mic(i) * 1.e-12
            gcdist    = gcdist    + scaling * cdisturb(i) * 1.e-12
            gtotfall  = gtotfall  + scaling * (falll(i) + fallr(i) + fallw(i)) * 1.e-12

         end do ! mlpt
      end do ! npoi1

      gtrunoff  = gsrunoff  + gdrainage
      gprcp     = gprcp     / gtarea
      gaet      = gaet      / gtarea
      gt        = gt        / gtarea
      gtrunoff  = gtrunoff  / gtarea
      gsrunoff  = gsrunoff  / gtarea
      gdrainage = gdrainage / gtarea
      gdwtot    = gdwtot    / gtarea
      aratio   = gaet     / gprcp
      rratio   = gtrunoff / gprcp
      sratio   = gsrunoff / max(dble(0.01),gtrunoff)
      tratio   = gt       / gaet

! write some diagnostic output to history file if I am the first proccess
      if (myid .eq. 0) then

         open(10,file=trim(outdir)//'/inland.log',status='unknown')

         if ( iyear0 .eq. iyear ) then
            write (10,*) ' '
            write (10,*) '**************************************'
            write (10,*) '* InLand - Surface model for the MBSCG'
            write (10,*) '**************************************'
            write (10,*) '* Version ',PACKAGE_VERSION,', November 2012'
#ifdef SINGLE_POINT_MODEL
            write (10,*) '* 0D: Single Point Model version.'
#endif /* SINGLE_POINT_MODEL */
            write (10,*) '**************************************'
            write (10,*) ' '

            write (10,*) 'Domain parameters ' 
            write (10,80) xres, yres
80          format (1x,'INFO: model lon, lat resolution (degrees) : ',2f8.2)
            write (10,90) nlon, nlat
90          format (1x,'INFO: model domain (nlon x nlat)          :     ',i3,'  x  ',i3) 
            write (10,95) snorth, ssouth, swest, seast
95          format (1x,'INFO: model limits (N,S,W,E)              : ',f8.2,' ',f8.2,' ',f8.2,' ',f8.2)
            write (10,*) ' '

            write (10,*) 'Temporal parameters '
            write (10,*) 'INFO: first year run in this sequence     : ', iyear0
            write (10,*) 'INFO: last year run in this sequence      : ', nrun + iyrlast
            write (10,*) 'INFO: length of this simulation (years)   : ', nrun
            write (10,*) 'INFO: number of iterations per day        : ', int(86400.0 / dtime)
            write (10,*) ' '

            write (10,*) 'Simulation parameters'
            write (10,*) 'INFO: Vegetation method                   : ',isimveg
            write (10,*) 'INFO: Fire method                         : ',isimfire
            write (10,*) 'INFO: CO2 method                          : ',isimco2
            write (10,60) co2init
60          format (1x,'INFO: initial co2 concentration in mol/mol: ',f8.6)
            write (10,70) o2init
70          format (1x,'INFO: initial o2 concentration in mol/mol : ',f8.6)
            write (10,*) 'INFO: Spin-up soil method                 : ',soilcspin
            write (10,*) 'INFO: Land Use method                     : ',isimland
            write (10,*) 'INFO: Infltration method                  : ',isinfilt
            write (10,*) 'INFO: Root water uptake method            : ',isimrwu
#ifndef SINGLE_POINT_MODEL
            write (10,*) ' '
            write (10,*) 'Heterogeneous Parameterization'
            write (10,*) 'INFO: Tauwood0 parameters                 : ',itauw
            write (10,*) 'INFO: Vmax parameters                     : ',ivmax
            write (10,*) 'INFO: Specla parameters                   : ',isla
            write (10,*) 'INFO: Carbon Allocation to wood parameters: ',ica
#endif

            if (isimagro .eq. 1) then
               write (10,*) ' '
               write (10,*) 'Crop Simulation parameters'
               write (10,*) 'INFO: Agro method                         : ',isimagro
               write (10,*) 'INFO: Crop type method                    : ',icroptype
            endif

            write (10,*) ' '
            write (10,*) 'I/O Parameters ' 
            write (10,*) 'INFO: input directory                     : ',indir
            write (10,*) 'INFO: input file                          : ',infile
            write (10,*) 'INFO: data directory                      : ',datadir
            write (10,*) 'INFO: output directory                    : ',outdir

         endif

         write (10,*) ' '
         write (10,*) '* * * annual diagnostic fields * * *',iyear
         write (10,*) ' '
         write (10,9001) gnee
         write (10,9000) gnpp
         write (10,9002) ggpp
         write (10,9010) gbiomass
         write (10,9020) galitc
         write (10,9021) gblitc
         write (10,9030) gsoic
         write (10,9032) gco2soi
         write (10,9034) galitn
         write (10,9036) gblitn
         write (10,9038) gsoin
         write (10,*) ' '
         write (10,9040) gprcp
         write (10,9050) gaet
         write (10,9060) gt
         write (10,9070) gtrunoff
         write (10,9080) gsrunoff
         write (10,9090) gdrainage
         write (10,9095) gdwtot
         write (10,*) ' '
         write (10,9100) aratio
         write (10,9110) rratio
         write (10,*) ' '
         write (10,9120) tratio
         write (10,9130) sratio
         write (10,*) ' '
!        call flush (200)
         write (201,9500) iyear, gnee, gnpp, gtotfall, gcdist, gco2mic,        &
                          gleach, ggpp, gbiomass, gsoic, gsoin, galitc,        &
                          gblitc, gco2soi, gprcp, gaet, gt, gsrunoff, gdrainage
!        call flush (201)

         write (STDOUT,*) ' '
         write (STDOUT,*) '* * * annual diagnostic fields * * *',iyear
         write (STDOUT,*) ' '
         write (STDOUT,9200) gtarea/1000.0/1000.0
         write (STDOUT,*) ' '
         write (STDOUT,9001) gnee
         write (STDOUT,9000) gnpp
         write (STDOUT,9002) ggpp
         write (STDOUT,9010) gbiomass
         write (STDOUT,9020) galitc
         write (STDOUT,9021) gblitc
         write (STDOUT,9030) gsoic
         write (STDOUT,9032) gco2soi
         write (STDOUT,9034) galitn
         write (STDOUT,9036) gblitn
         write (STDOUT,9038) gsoin
         write (STDOUT,*) ' '
         write (STDOUT,9040) gprcp
         write (STDOUT,9050) gaet
         write (STDOUT,9060) gt
         write (STDOUT,9070) gtrunoff
         write (STDOUT,9080) gsrunoff
         write (STDOUT,9090) gdrainage
         write (STDOUT,9095) gdwtot
         write (STDOUT,*) ' '
         write (STDOUT,9100) aratio
         write (STDOUT,9110) rratio
         write (STDOUT,*) ' '
         write (STDOUT,9120) tratio
         write (STDOUT,9130) sratio
         write (STDOUT,*) ' '
      endif
 
9200  format (1x,'total area            of the domain (km^2)    : ', f24.1)
9000  format (1x,'total npp             of the domain (gt-c/yr) : ', f24.14)
9001  format (1x,'total nee             of the domain (gt-c/yr) : ', f24.14)
9002  format (1x,'total gpp             of the domain (gt-c/yr) : ', f24.14)
9010  format (1x,'total biomass         of the domain (gt-c)    : ', f24.14)
9020  format (1x,'aboveground litter    of the domain (gt-c)    : ', f24.14)
9021  format (1x,'belowground litter    of the domain (gt-c)    : ', f24.14)
9030  format (1x,'total soil carbon     of the domain (gt-c)    : ', f24.14)
9032  format (1x,'total soil co2 flux   of the domain (gt-c)    : ', f24.14)
9034  format (1x,'aboveground litter n  of the domain (gt-c)    : ', f24.14)
9036  format (1x,'belowground litter n  of the domain (gt-c)    : ', f24.14)
9038  format (1x,'total soil nitrogen   of the domain (gt-c)    : ', f24.14)
9040  format (1x,'average precipitation of the domain (mm/yr)   : ', f24.14)
9050  format (1x,'average aet           of the domain (mm/yr)   : ', f24.14)
9060  format (1x,'average transpiration of the domain (mm/yr)   : ', f24.14)
9070  format (1x,'average runoff        of the domain (mm/yr)   : ', f24.14)
9080  format (1x,'average surf runoff   of the domain (mm/yr)   : ', f24.14)
9090  format (1x,'average drainage      of the domain (mm/yr)   : ', f24.14)
9095  format (1x,'average moisture recharge of domain (mm/yr)   : ', f24.14)
9100  format (1x,'total aet      / precipitation                : ', f24.14)
9110  format (1x,'total runoff   / precipitation                : ', f24.14)
9120  format (1x,'transpiration  / total aet                    : ', f24.14)
9130  format (1x,'surface runoff / total runoff                 : ', f24.14)
9500  format (1x,i4,18f24.14)

      return
end subroutine gdiag 
