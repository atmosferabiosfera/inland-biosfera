#include "inland_config.h"
!
! #    #    ##       #    #    #
! ##  ##   #  #      #    ##   #
! # ## #  #    #     #    # #  #
! #    #  ######     #    #  # #
! #    #  #    #     #    #   ##
! #    #  #    #     #    #    # OFFLINE!
!
! ---------------------------------------------------------------

program main !main_offline!
! ---------------------------------------------------------------
! common modules
#ifndef SINGLE_POINT_MODEL
      use inland_parameters
      use inland_control, only: jday, iday, imonth, iyear, iyear0, idailyout, &
                                imonthout, iyearout, isimveg, isimfire,       &
                                isimco2, nspinsoil, spinmax, eqyears, istep,  &
                                spincons, spinfrac, ipointout, ccmexist, & 
                                nluc, iyrluc, isimland, &
                                indir, infile, datadir, outdir, &
                                env_ran2val, env_fastexec, env_debug, &
                                vegtypefile, hrmapfile, itauw,isla,ica,ivmax, &
                                isinfilt, isimrwu, cropsfile, nrun
      use inland_comatm, only: co2conc, o2conc
      use inland_subgrid
#else /* SINGLE_POINT_MODEL */
!     TODO: verify if there are variables that should be no longer used - fzm
      use inland_parameters, only: npoi, nsoilay, ndaypm, dtime, hvap, lbeg,   &
                                 lend, mpt, ndaypy, nlon, nlat, rhow, xres,  &
                                 yres
      use inland_control, only: jday, iday, imonth, iyear, iday0, imonth0,   &
                                iyear0, idailyout, imonthout, iyearout,      &
                                isimveg, isimfire, isimco2, nspinsoil, istep,&
                                spinmax, eqyears, spincons, spinfrac,        &
                                ipointout, ccmexist, nan, &
                                indir, infile, datadir, outdir, &
                                env_ran2val, env_fastexec, env_debug, &
                                vegtypefile, hrmapfile, isinfilt, isimrwu, &
                                indir, infile, datadir, outdir, nrun, cropsfile
      use inland_comforc, only: xprec, xcld, xlati, xlin, xsin, xta, xqa,    &
                                xua, dimforc
#endif /* SINGLE_POINT_MODEL */
      use inland_com1d
      use inland_comatm
      use inland_comcrop
      use inland_comhour
      use inland_comdiag
      use inland_comdiag, only: diagstart, diagend, ndiagpt, nfreq
      use inland_comwork
      use inland_comfire
      use inland_lsmmapib, only: lati, loni
      use inland_comveg, only: exist, topparu, topparl, tneetot,plai,tl, use, plaimx, aleaf, aroot, awood, cil4, fallrsgc
      use inland_combgc, only:cnroot, cnleaf,rconst, cnwood, cnroot
      use inland_comhyd, only:gsuvap,gtrans
      use inland_comsoi, only: tsoi, soihfl, h20, smsoil, cndepth, wsoi, fwpud

      implicit none

! local variables
      integer istyrd,      & ! first yr daily means exist (from precip file)
              !istyrm,     & ! first yr monthly anoms exist (from file)
              iwest,       & ! 1st lon index for subset
              iy1,         & ! first year for year loop
              iy2,         & ! last year for year loop
              iyrdaily,    & ! year to start reading daily means
              iyrlast,     & ! last year of previous run (for restart)
              iyrmon,      & ! month year of previous run (for restart)
! The below variable has been inherited from CCSM3-ibis.
              imonthlast,  & ! last month of previous run (for restart)
              idiag,       & ! number of diagnostic files requested
              irestart,    & ! 0: normal mode 1: restart mode
              iholdsoiln,  & ! 0: don't save inorganic soil N values 1: keep inorganic soil N values 
              iyrrestart,  & ! Year to restart model
              jnorth,      & ! 1st lat index for subset
              nanom,       & ! # of years in the anomaly files
              nanomd,      & ! # of years in the daily anomaly files
              niter,       & ! total number of time iterations per day
              nday,        & ! number of days since the start of the simulation
              plen,        & ! length of precipitation event in timesteps (see plens)
              plenmax,     & ! upper bound on plen
              plenmin,     & ! lower bound on plen

              seed,        & ! first value used to initialize the random # generator
              seed2,       & ! second value used to initialize the random # generator
              seed3(32),   & ! third value used to initialize the random # generator
              seed4,       & ! fourth value used to initialize the random # generator

              spin,        & ! counter for iterations used to spin up soilbgc
              soilcspin,   & ! 0: no spinup procedure for soil c  1: acceleration procedure used
              irrigate,    & ! 0: no irrigation 1: irrigation used
              lun, ifile,  & ! file indices
              timed,       & ! time of day since midnight (in seconds)
              i, j, n, ivar   ! loop indices

      integer daystart       ! index of when to start the day loop
      integer daystop        ! month to end simulation (12 by default)
      integer monthstop      ! month to end simulation (12 by default)

! Imbuzeiro: linenum is created to print the input files (xta, xqa, xua, xprec,
!            xsin, xlin) as an output
      integer istat

      real *8 rn,          & ! net radiation flux (short and long waves)
              pari,        & ! incoming PAR
              apar,        & ! APAR
              paro,        & ! outgoing PAR
              reflsw,      & ! reflected solar radiation
              swin,        & ! incoming solar radiation (W/m²)
              swout          ! reflect solar radiation (W/m²)

#ifdef SINGLE_POINT_MODEL

! Imbuzeiro: linenum is created to print the input files (xta, xqa, xua, xprec,
!            xsin, xlin) as an output
      integer linenum,     & ! line number
              monthstart,  & ! index of when to start the month loop
              isoilforc,   & ! 0: fixed soil, 1: dynamic soil
              timel,       & ! local time of day since midnight (in seconds)
              nc             ! 0D specific loop indice

      real*8  soilmoistl
! Hewlley's variables as of june, 2007
! FIXME: this is no place for this kind of variables.
!       -fzm
!      real*8  swnet,       & ! net shortware radiation
!              lwnet,       & ! net longwave radiation
!              vegt,        & ! vegetation canopy temperature
!              avgsurft,    & ! average surface temperature
!              surfalb,     & ! surface albedo
!              soilwet,     & ! total soil wetness
!              canopint,    & ! total canopy water storage
!              rootmoist,   & ! root zone soil moisture
!              canheat,     & ! canopy heat storage (J.m^-2)
!              canheatl,    & ! last canopy heat storage (J.m^-2)
!              delcanheat,  & ! change incanopy heat storage (J.m^-2)
!!              delsoilmoist,& ! change in soil moisture (kg.h2o.m^-2)
!              soilmoistsum,& ! soil moisture in all layers (kg.h2o.m^-2)
!              soilmoist(nsoilay), & ! soil moisture in which layer (kg.h2o.m^-2)
!              soilmoistl,  & ! last soil moisture (kg.h2o.m^-2)
!              autoresp,    & ! autotrophic respiration (kg.m^-2.s^-2)
!              surfheat,    & ! surface heat storage (J.m^-2)
!              delsurfheat, & ! change in surface heat storage (J.m^-2)
!              surfheatl,   & ! last surface heat storage (J.m^-2)
!!              delsurfstor, & ! change in surface heat storage (J.m^-2)
 !             wtotl,    & ! total wather stored in soil, vegetation and snow
!                            ! at previous timestep (kg.h2o.m^-2.s^-1)
 !             totsoilcarb, & ! total soil carbon (kg.C.m^-2)
! FIXME: surfalb or albedo? are them the same? -fzm
  !            albedo,      & ! surface albedo
  !            solartot,    & ! total incoming radiation (dir+diff,vis+nIR)
! FIXME: avgsurft or soiltemp? are them the same? -fzm
    !          soiltemp,    & ! average layer soil temperature (K)
   !           undefined,   & ! undefined by albedo
         !     delwtot,     & ! change in total water stored (kg.h2o.m^-2.s^-1)
     !         fpar,        & ! absorbed fraction of PAR
      !        totlai,      & ! total ecosystem lai
       !       totlivbio,   & ! total living biomass (kg.C.m^-2)
        !      agb            ! above ground biomass
#endif /* SINGLE_POINT_MODEL */

      real*8  co2init,     & ! initial co2 concentration in mol/mol
              o2init,      & ! initial o2 concentration in mol/mol
              dran,        & ! random number from generator
              plens,       & ! length of precipitation event in seconds
              startp,      & ! time to start precipitation event (seconds since midnight)
              endp,        & ! time to end precipitation event (seconds since midnight)
              ilens,       & ! length of irrigation event in seconds
              starti,      & ! time to start irrigation event (seconds since midnight)
              endi,        & ! time to end irrigation event (seconds since midnight) 
              slope,       & ! rate of decrease of number of iteration for soil C spinup 
              snorth,      & ! north latitude for subsetting std grid
              ssouth,      & ! south latitude for subsetting std grid
              swest,       & ! west longitude for subsetting std grid
              seast,       & ! east longitude for subsetting std grid
              ffact,       & ! numerical multiplying factor applied to N fertilizer being applied after 2000 
              test,        & ! test on timestep
              thigh          ! tower elevation (inland.infile) 
	      
      character*255  domain
      character*1024 filen
      character*4    chyear
      character*255  chtmp

#ifndef SINGLE_POINT_MODEL
! External
      real ran2              ! Function : random number generator
#endif /* SINGLE_POINT_MODEL */

! default values for subset grid
      data snorth, ssouth, swest, seast / 90., -90., -180., 180. /

! Ibis - CCSM3 variables
      logical rdlsf        ! set initial data for Antarctica, Greenland if true
      integer loopi,      &! index of little vector in big vector
              kpti,       &! start of vector segment
              kptj         ! end of vector segment

! namelist declaration
#ifndef SINGLE_POINT_MODEL
      namelist/INLAND_GRID/ irestart, iyrrestart, iyear0, nrun, iyrdaily, iyrmon, dtime, soilcspin, isimveg, isimfire, isimco2, co2init, o2init, isinfilt, isimrwu, isimland, iyrluc, nluc, mlpt, vegtypefile, hrmapfile, cropsfile, iyearout, imonthout, idailyout, idiag, imetyear, dmetyear, imetend, dmetend, irrigate, isimagro, icroptype, iwheattype, cropsfile, irotation, iholdsoiln, ffact, isoilay, elevin, thigh, domain, itauw, ivmax, isla,ica, rddaydims, isparse
#else
      namelist/INLAND_SINGLE_POINT/ iyear0, nrun, dtime, soilcspin, isimveg, isimfire, isoilforc, isimco2, co2init, o2init, isinfilt, isimrwu, snorth, ssouth, swest, seast
#endif /* SINGLE_POINT_MODEL */

! multiproc variables
      integer ithread, numthreads 
#ifdef _OPENMP
#ifndef SINGLE_POINT_MODEL
      !     USE omp_lib
      INTEGER OMP_GET_MAX_THREADS
      INTEGER OMP_GET_NUM_THREADS
      INTEGER OMP_GET_NUM_PROCS
      INTEGER OMP_GET_THREAD_NUM     
#endif
#endif

! time report variables
#ifdef __GFORTRAN__
      real, dimension(2) :: tarray
      real tresult
      integer telapsed
      telapsed = time()
#endif

      ipointout=0

! ---------------------------------------------------------------
!                       j o b   c o n t r o l 
! ---------------------------------------------------------------
!
! multiproc variables - put default values for now
      ithread = 1
      loopi = 1
      numthreads = 1
      numlv = 1

! model I/O variables
      call rd_env

! open input file (namelist)
! TODO: add option to specify infile file name in command line.
! FIXME: Verify if every variable is allocated or it is needed another
!      'inland_alloc' call here.

! setup defaults - document any changes here in inland-grid.namelist and inland-single_point.namelist
! temporal parameters
      irestart   = 0        ! 0: not a restart run  1: restart run (default 0)
      iyrrestart = 9999     ! year to restart model (default 9999)
      iyear0     = 1981     ! initial year of simulation (don't change for restart) (default 1981)
      nrun       = 2        ! number of years in this simulation (change for restart) (default 2)
      iyrdaily   = 9999     ! year to start reading daily data (ditto) (default 9999)
      iyrmon     = 9999     ! year to start reading montly data (default 9999)
      dtime      = 3600     ! time step in seconds (default 3600)
!
! simulation parameters
      soilcspin  = 0        ! 0: no soil spinup, 1: acceleration procedure used (default 0)
      isimveg    = 0        ! 0: static veg, 1: dynamic veg, 2: dynamic veg-cold start (default 0)
      isimfire   = 1        ! 0: no fire disturbance (0%/yr), 1: natural const (0.5%/yr), 2: CTEM, 3: IBIS (default 1, ignored if isimveg=0)
      isimco2    = 0        ! 0: fixed co2,  1: ramped co2 (default 0)
      co2init    = 0.000350 ! initial co2 concentration in mol/mol (real) (default 0.000350)
      o2init     = 0.209000 ! initial o2 concentration in mol/mol (real) (default 0.209000)
      isinfilt   = 0        ! isinfilt - Infiltration Function - 0: according to Darcy (1856); 1: according to Green-Ampt (1911) (default 0)
      isimrwu    = 0        ! isimrwu - Root water uptake module - 0: according to Foley et al., 1996; 1 according to Li et al. (2006) (default 0)
      mlpt       = 1        ! multiplicity of land points (subgrid tiles) (default 1, no tiles)
#ifndef SINGLE_POINT_MODEL
      isimland   = 0        ! 0: fixed land use,  1: dynamic land use
      iyrluc     = 999      ! year to start reading land transition data
      nluc       = 50       ! number of years to read land transition data
      vegtypefile= ''       ! file to read initial vegtype (and tileprop if mlpt>1) (default '', read default)
      hrmapfile  = ''       ! read high-res data to produce initial veg map (default '', don't read)
#endif
!TODO put variables for agro in separate namelist
      cropsfile  = ''       ! read crops data to produce crop veg map (default '', don't read)
!

! output parameters
      iyearout   = 1        ! 0: no yearly output, 1: yearly output (default 1)
      imonthout  = 1        ! 0: no monthly output, 1: monthly output (default 1)
      idailyout  = 0        ! 0: no daily output, 1: daily output (default 0)
      idiag      = 0        ! 0: no diagnostic output, 1-10 # of files to output (default 0)
!
! crops simulation paramters
imetyear   = 9999      ! year to start using hourly met data
dmetyear   = 1         ! day to start using hourly met data (17 geovane)
imetend    = 9999      ! year to end reading in hourly met data
dmetend    = 119       ! day to end reading in hourly met data (39 geovane)
irrigate   = 0         ! 0: no irrigation 1: irrigation strategy used
isimagro   = 0         !  isimagro - 0: no agro crops, 1: unique crop defined by icroptype, 2: crops defined by data in cropsfile
icroptype  = 0         ! define crop type to use in all grid points - 13: soybean / 14: maize / 15: wheat / 16: sgc
iwheattype = 1         ! 1: spring wheat 2: winter wheat - used if icroptype = wheat
cropsfile  = ''        ! cropsfile - read crop data to produce crop veg map (use mlpt>1)
!
irotation  = 0         ! 0: none 1: winter wheat/fallow 2: 2 corn/soy 3: corn/soy/spring wheat 4: soy/winter wheat/corn
iholdsoiln = 1         ! 0: doesn't save soil inorganic N from restart 1: save inorganic soil N  
ffact      = 1.0       ! numeric multiplying factor applied to N fertilizer after 2000 (for all crops)
isoilay    = nsoilay-2 ! soil layer for which nitrate leaching/drainage is output  
elevin     = 550       ! site elevation        (m)
thigh      = 6         ! tower input elevation (m)
!
! domain parameters
domain     = ''
snorth     =  -18.75   ! snorth - northern latitude for subsetting in/output (no default)
ssouth     =  -25.75   ! ssouth - southern latitude for subsetting in/output (no default)
swest      =  -53.75   ! swest - western longitude for subsetting in/output (no default)
seast      =  -43.75   ! seast - eastern longitude for subsetting in/output (no default)
!
#ifndef SINGLE_POINT_MODEL 
!
! Heterogeneous Parameterization
      itauw	 =  0	! itauw -  0:reads 1 dimension tauwood0 from parms canopy; 1:reads 2 dimension from input maps
      ivmax 	 =  0   ! ivmax -  0:reads 1 dimension vmax_pft from parms canopy; 1:reads 2 dimension from input maps
      isla       =  0   ! isla  -  0:reads 1 dimension specla from parms canopy; 1:reads 2 dimension from input maps
      ica        =  0   ! ica   -  0:reads 1 dimension carbon allocation to wood, leaf and root from parms canopy; 1:reads 2 dimension from input maps
!
#endif

! different defaults for single_point
#ifdef SINGLE_POINT_MODEL
      isoilforc  = 0        ! 0: dynamic soil  physics 1: forced soil physics (default 0)
      iyear0     = 1992
      nrun       = 12
      iyearout   = 0
      imonth     = 0
      isimveg    = 1
#endif

!read namelist
       lun = 12
      open (lun, status='old', file=trim(infile))
#ifndef SINGLE_POINT_MODEL
      read (lun, end=99, NML=INLAND_GRID)
#else
      read (lun, end=99, NML=INLAND_SINGLE_POINT)
#endif
99    close(lun)

! if domain name was given, look it up in conf/domains file
      if ( trim(domain) .ne. '' ) then

         filen=trim(indir)//'/'//'conf/inland-grid.domains'
         open(15, status='old', file=filen, err=8102)

8001     read(15,*,end=8152,err=8103) chtmp,snorth,ssouth,swest,seast
         ! got to end, raise error
         if (trim(chtmp) .eq. 'stop') then
            goto 8151
         endif
         ! found domain name, exit loop
         if (trim(chtmp) .eq. trim(domain)) then
            goto 8152
         endif
         goto 8001

         ! error handling
8103     call commenthandler (15, filen)
         goto 8001
8102     write(*,*) 'ERROR: ibis.infile domain='//trim(domain)//' but conf/inland-grid.domains file not found'
         stop 1
8151     write(*,*) 'ERROR: ibis.infile domain='//trim(domain)//' but conf/inland-grid.domains does not contain '//trim(domain)
         stop 1
       
8152     close (15)

      end if


! why was this in the inland.infile section?
#ifdef SINGLE_POINT_MODEL
      linenum = 0
#endif

! we are not coupled. ccmexist (from inland_control.F90) is 0.
      ccmexist = 0
      
! Check values in inland-grid.infile
#ifdef SINGLE_POINT_MODEL
      call check(soilcspin, nrun, isoilforc)
#else
      call check(irestart, soilcspin, nrun, snorth, ssouth, swest, seast)      
#endif

! tell user about the simulation
      write (*,*) ' '
      write (*,*) '**************************************'
      write (*,*) '* InLand - Surface model for the MBSCG'
      write (*,*) '**************************************'
      write (*,*) '* Version ',PACKAGE_VERSION,', November 2012'

#ifdef SINGLE_POINT_MODEL
      write (*,*) '* 0D: Single Point Model version.'
#endif /* SINGLE_POINT_MODEL */
      write (*,*) '**************************************'
      write (*,*) ' '
      if (irestart .eq. 1) then
         write (*,*) 'INFO: running in restart mode'
#ifdef SINGLE_POINT_MODEL
         write (*,*) 'ERROR: Single Point Model (0D) does NOT support restart.' 
         write (*,*) 'Please turn off restart (set irestart to 0) in config file.'
         stop 1
#endif /* SINGLE_POINT_MODEL */
      end if

      write (*,*) ' '
      write (*,9010) xres, yres
      write (*,9020) nlon, nlat
      write (*,9025) snorth, ssouth, swest, seast
      write (*,*) ' '
      write (*,9041) trim(indir)
      write (*,9042) trim(infile)
      write (*,9043) trim(datadir)
      if ( trim(outdir) .ne. '.' ) then
         write (*,9044) trim(outdir)
      end if
#ifndef SINGLE_POINT_MODEL
      if ( env_ran2val .ne. 0. ) then
         write (*,9045) env_ran2val
      end if
#endif /* SINGLE_POINT_MODEL */
      write (*,*) ' '
      write (*,9000) nrun

9000  format (1x,'INFO: length of this simulation (years)   : ',i8)
9001  format (1x,'INFO: Vegetation method                   : ',i5)
9002  format (1x,'INFO: CO2 method                          : ',i5)
9006  format (1x,'INFO: number of iterations per day        : ',i8)
9010  format (1x,'INFO: model lon, lat resolution (degrees) : ',2f8.2)
9020  format (1x,'INFO: model domain (nlon x nlat)          :     ',i3,'  x  ',i3)  
9025  format (1x,'INFO: model limits (N,S,W,E)              : ',f8.2,' ',f8.2,' ',f8.2,' ',f8.2)  
9030  format (1x,'INFO: last year run in this sequence      : ',i8)
9040  format (1x,'INFO: number of threads                   : ',i8)
9041  format (1x,'INFO: input directory                     :     ',a)
9042  format (1x,'INFO: input file                          :     ',a)
9043  format (1x,'INFO: data directory                      :     ',a)
9044  format (1x,'INFO: output directory                    :     ',a)
9045  format (1x,'INFO: randomval                           :     ',f7.4)
9055  format (1x,'NOTICE - setting numthreads to ',i2,' instead of OMP_GET_MAX_THREADS ( ',i2,' )')

      if (iyrdaily.eq.iyrmon .and. iyrdaily.ne.9999) then
          write(*,*) 'ERROR: cannot run because iyrdaily:',iyrdaily,' equal irymon:', iyrmon
          stop 2
      endif

! ---------------------------------------------------------------
!                     t i m e   c o n t r o l 
! ---------------------------------------------------------------
! determine the number of timesteps per day
      niter = int(86400.0 / dtime)

! test on length of the time step and total number of iteration per day
      test = dmod(dble(86400.0), dtime)

      if (test .gt. 1.e-20) then
         write (*,*) 'ERROR: dtime ', dtime, ' should be divisible into 86400'  
         stop 3
      else
         write (*,9006) niter
      end if

! lbeg, lend, mpt, kti and kptj now set in readit subroutine
! ---------------------------------------------------------------
!                   i n i t i a l i z a t i o n
! ---------------------------------------------------------------
! Allocate variables needed by rd_param and readit (except those with npoi dim).
      call inland_prealloc

! check for restart conditions and set "iyrlast" accordingly, open
! inland in append mode
      if (irestart .eq. 1) then
         write(chyear(1:4),'(i4)') iyrrestart
         filen = trim(outdir)//'/yearsrun.'//chyear//'.dat'
         open(13, status='old', file=filen, err=9050)
         read(13,*) iyrlast
! those variables to turn possible restart of ran2 subroutine
         read(13,*) seed
         read(13,*) seed2
         read(13,*) seed3
         read(13,*) seed4
         read(13,*) aleafi
         read(13,*) astem
         read(13,*) astemi
         read(13,*) arepr
         read(13,*) iday
         read(13,*) cnleaf
         read(13,*) cnwood
         read(13,*) cnroot
         read(13,*) h20
         read(13,*) cndepth
         read(13,*) vf
         read(13,*) iniday
         read(13,*) fleaf
         read(13,*) fallrsgc
         close(13)
         goto 9059
  9050   write(*,*) 'ERROR: inland-grid.namelist restart flag=1, indicating restart, but output/yearsrun.',iyrrestart,'.dat not found'
         stop 4
         write(*,*) 'beginning from ', iyear0
         iyrlast = iyear0 - 1
         open(20,file='inland.out.global',status='unknown',access='append')
         open(30,file='inland.out.vegtype',status='unknown',access='append')
      else
         iyrlast = iyear0 - 1
9059  end if

      write (*,9030) nrun + iyrlast
      write (*,*) ' '

! call RD_PARAM to read in all parameter files
      call rd_param(irestart)

! initial concentration of co2, o2
      co2conc = co2init
      o2conc  = o2init

! read global boundary condition datasets
      call readit(isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
#ifndef SINGLE_POINT_MODEL
      if(npoi .eq. 1) then
         call build_file
      endif
#endif /* SINGLE_POINT_MODEL */
      
! multiproc variables
#ifdef _OPENMP
#ifndef SINGLE_POINT_MODEL
      numthreads = OMP_GET_MAX_THREADS() 

      ! make sure numthreads <= 2, unless explicitly set with OMP_NUM_THREADS
      call getenv("OMP_NUM_THREADS", chtmp)
      if ( numthreads .gt. 2 .and. trim(chtmp) .eq. "" ) then
         print *, ''
         write (*,9055) 2, OMP_GET_MAX_THREADS()
         numthreads = 2
         numlv = numthreads
         call OMP_SET_NUM_THREADS(numthreads)
      end if

      ! make sure numthreads <= npoi
      if ( numthreads > npoi ) then
         print *, ''
         write (*,9055) npoi, OMP_GET_MAX_THREADS()
         numthreads = npoi
         call omp_set_num_threads(numthreads)
      end if
      numlv = numthreads
      if ( numthreads .gt. 1 ) write (*,9040), numthreads 
#endif
#endif

      call alloc(irestart)

#ifdef SINGLE_POINT_MODEL
! read forcing data
      call readforc(isoilforc)

! FIXME: 'test' variable. No clue what it is for, this is a badly written code.
      test = dimforc
      
      call build_file      

#endif /* SINGLE_POINT_MODEL */

! map coupled model variables based on variables read from data files
      call map_offline(imonthlast, rdlsf, loopi, kpti, kptj, ipointout)

! check if diagnostic output is requested, if so read info from 'diag.infile'
! TODO: not supported yet, can be added as a later feature. Just keep '0' at
!       conf/inland-<model>.infile for 'idiag'.
      !if (idiag .ne. 0) call inidiag(idiag)

#ifndef SINGLE_POINT_MODEL
! preliminary analysis of climate data
      if (irestart.eq.0 .or. isimveg.eq.0) then
         call climanl
      end if

      if (iyrdaily .le. iyrlast+nrun) then
         write(chyear(1:4),'(i4)') iyrdaily
         filen = trim(datadir)//'/daily/prec.daily.'//chyear//'.nc'
         call inird(filen,istyrd)
         if (istyrd .le. 0) istyrd = 9999
      end if

!    If first year of this run/restart is an anomaly year, then read in 
! anomalies for month 12 of previous year and month 1 of this year
! (note: rdanom reads imonth+1) 
      iy2 = iyrlast + nrun
#endif /* SINGLE_POINT_MODEL */

! initialize the model
      call initialib(isimveg, irestart, iyrlast, imonthlast, rdlsf, lati, loni)
#ifndef SINGLE_POINT_MODEL
! initialize crop variables
! have to be initialized when crops are replacing natural vegetation in
! a series of simulations
      if (isimagro .gt. 0) call initialcrop(irestart)

! initialize random number generator seed
      if(irestart.eq.0) then
         seed = -1
         seed2 = 123456789
         seed3 = 0
         seed4 = 0
         dran = ran2 (seed,seed2,seed3,seed4)
      endif

! make sure all points have a unique seed
!      do i = lbeg, lend
!         seedvec(i) = i
!      end do
! make sure subgrids get same seed at given point, this to have equal random values 
      do i = 1,npoi1
         do j = 1, mlpt
            seedvec( subgrid_get_index(i,j) ) = i
         end do
      end do

#endif /* SINGLE_POINT_MODEL */

! ---------------------------------------------------------------
!              s t a r t   o f   t i m e   l o o p s
! ---------------------------------------------------------------
      write (*,*) ' '
      write (*,*) '********************************'
      write (*,*) '* start of the INLAND simulation *'
      write (*,*) '********************************'
      write (*,*) ' '
! reset elapsed number of days, accumulating past days if restart mode
      nday = 0

      if (irestart .eq. 1) then
         do 150 iyear = iyear0, iyrlast
            nday = nday + 365
            if (mod(iyear,4).eq.0) then
               if (mod(iyear,100).ne.0) then
                  nday = nday + 1
               elseif (mod(iyear/100,4).eq.0) then
                  nday = nday + 1
               end if
            end if
150      continue
      end if

!---------------------- YEAR YEAR YEAR -------------------------------------
! start of yearly loop
!---------------------- YEAR YEAR YEAR -------------------------------------
      iy1 = iyrlast + 1
      iy2 = iyrlast + nrun
      do 200 iyear = iy1, iy2
! reset julian date
         jday = 0

!    Determine the calendar for this year.
!    Leap years occur every calendar year that is divisible by 4 (such as
! 1980, 1984, etc.) *except* for years which begin centuries not divisible
! by 400
!    For example, 1700, 1800, and 1900 are *not* leap years, but 1600
! and 2000 *are* leap years
         ndaypm(2) = 28
         ndaypy = 365

! 0D INLAND does not consider leap years.
#ifndef SINGLE_POINT_MODEL
         if (mod(iyear,4).eq.0) then
            if (mod(iyear,100).ne.0) then
               ndaypm(2) = 29
               ndaypy = 366
            elseif (mod(iyear/100,4).eq.0) then
               ndaypm(2) = 29
               ndaypy = 366
            end if
         end if
#endif /* SINGLE_POINT_MODEL */
! NOTICE: Although Hewlley's version has CO2 commented out, you can disable it
!        just by not setting isimco2=1 on the infile.
! get new o2 and co2 concentrations for this year
         if (isimco2.eq.1) call co2(co2init)

!---------------------- MONTH MONTH MONTH ----------------------------------
! start of monthly loop
!---------------------- MONTH MONTH MONTH ----------------------------------
         if ( env_fastexec .gt. 0 ) then
            monthstop = 1
         else
            monthstop = 12
         end if
#ifdef SINGLE_POINT_MODEL
! Observed, step data, might not start on first month of year.
         if (iyear.eq.iy1) then 
            monthstart = imonth0
         else
            monthstart = 1
         endif
         do 210 imonth = monthstart, monthstop
! rdanom is not called on single point model.
#else /* SINGLE_POINT_MODEL */
         do 210 imonth = 1, monthstop
             if (iyear.ge.iyrmon) then
     	        call rdmon(iyear, iwest, jnorth)
	     endif
#endif /* SINGLE_POINT_MODEL */

!------------------------- DAY DAY DAY -------------------------------------
! start of daily loop
!------------------------- DAY DAY DAY -------------------------------------
               daystart = 1
#ifdef SINGLE_POINT_MODEL
! Observed, step data, might not start on first day of the month for the first
!   month.
            if ( (iyear.eq.iy1).and.(imonth.eq.imonth0) ) then
               daystart = iday0
            endif
#endif /* SINGLE_POINT_MODEL */

         if ( env_fastexec .gt. 0 ) then
            daystop = 1
         else
            daystop = ndaypm(imonth)
         end if
!            do 220 iday = daystart, ndaypm(imonth)
            do 220 iday = daystart, daystop

! update user on model progress once a day
               write (STDOUT,9100,advance='no') char(13), char(13), iyear, imonth, iday
#ifdef IFORT_NOFLUSH_BUG
               close(STDOUT)
               open(STDOUT)
#endif
#ifdef PGF90_NOFLUSH_BUG
               call flush(STDOUT)
#endif
9100           format (a,79x,a,'Running now: ',i4.4,'/',i2.2,'/',i2.2)
9101           format (a,79x,a,'Running now: ',i4.4,'/',i2.2,'/',i2.2,a)
9102           format (a,79x,a,'Completed ',i4.4,'/',i2.2,'/',i2.2)

! update elapsed number of days and julian date
               nday = nday + 1
               jday = jday + 1

#ifdef SINGLE_POINT_MODEL
               call daily(test,dimforc)
#else /* SINGLE_POINT_MODEL */
    if(isimagro .gt. 0) then
! update crop calendar
        if(irestart .eq. 0) then
            do i = lbeg, lend
               do j = scpft, ecpft
                  if (exist(i,j).eq.1.and.iday.eq.pcd(j).and.imonth.eq.pcm(j)) then
                           ncyears(i)=ncyears(i)+1
                           cdays(i)=0
                           pstart(i,j)=999
                  endif
               enddo
            enddo
!
            do i = lbeg, lend
               if (cdays(i).lt.366) then
                  cdays(i) = cdays(i) + 1
               endif

            enddo
        else !if irestart
!#### TODO optimization this if!!!!
            if(iyear .eq. iyrrestart+1 .and. imonth .eq. 1 .and. iday .eq. 1)then
               do i = lbeg, lend
                 do j = scpft, ecpft
                    if (exist(i,j).eq.1.and.iday.eq.pcd(j).and.imonth.eq.pcm(j)) then
                             ncyears(i)=ncyears(i)+1
                             cdays(i) = cdays(i) + 1
                             pstart(i,j)=999
                    endif
                 enddo
               enddo
                 do i = lbeg, lend
                      if (cdays(i).lt.366) then
                         cdays(i) = cdays(i) + 1
                      endif
                 enddo
              else
               do i = lbeg, lend
                 do j = scpft, ecpft
                    if (exist(i,j).eq.1.and.iday.eq.pcd(j).and.imonth.eq.pcm(j)) then
                             ncyears(i)=ncyears(i)+1
                             cdays(i)=0
                             pstart(i,j)=999
                    endif
                 enddo
               enddo
               do i = lbeg, lend
                  if (cdays(i).lt.366) then
                     cdays(i) = cdays(i) + 1
                  endif
               enddo


            endif !end if iyrestart
        endif! end if irestart
     endif
! get daily means (if requested)
! and determine today's climatic conditions

! TODO: port the weather generator!
!      maybe we need to sed the kpti-kptj interval to allow parallelizing the
!      weather generator. by now, we just use lbeg:lend but keep your eyes open!
               if (iyear.ge.iyrdaily) then
                  call rdday(jday,istyrd,iwest,jnorth)
                  call daily(seed,seed2,seed3,seed4,1,iyrlast,nrun)
               else
                  if (iyear.ge.iyrmon) then
                     call dailymon(seed,seed2,seed3,seed4,0,iyrmon,iyrlast,nrun)
                  else
                     call daily(seed,seed2,seed3,seed4,0,iyrlast,nrun)
                  endif
               end if

#endif /* SINGLE_POINT_MODEL */

! determine the daily vegetation cover characteristics
               call pheno(kpti, kptj)


      if (isimagro .gt. 0) then

          if (irotation .gt. 0 .and. jday .eq. 1) call rotation(irestart, iyrrestart)

         call planting(irestart,iyrrestart,jday,ffact)
!
! call daily crop phenology subroutine
!

         call phenocrop(kpti,kptj)
      endif
!
! call soil biogeochemistry model
!
! soil spin up procedure calls soil bgc model spincons times
! for each time step to spin up the process (accelerate C & N pools)
               spinfrac  = 0.75     ! fraction of nspinsoil time used to
!                                   ! spin up soil at maximum spin up rate
!
               spincons  = 40.0     ! number of times soilbgc subroutine is
!                                   ! called for each day in simulation during
!                                   ! acceleration procedure
!
               eqyears   = 50       ! years used past spin up to slowly bring
!                                   ! entire soil to a steady equilibrium
!
               nspinsoil=iyear0+150 ! year of simulation that soil c reaches equilibrium

               if (soilcspin .eq. 1) then
                  if ((iyear-iyear0) .le. (spinfrac*(nspinsoil-iyear0-eqyears))) then
                     spinmax = int(spincons)
                  elseif ((iyear-iyear0) .lt. (nspinsoil-iyear0-eqyears)) then
                     slope   = spincons / ((nspinsoil - iyear0 - eqyears) -   &
                               (spinfrac * (nspinsoil - iyear0 - eqyears)))
                     spinmax = int(spincons - (slope * ((iyear - iyear0) -    &
                               (spinfrac * (nspinsoil - iyear0 - eqyears)))))
                     spinmax = max(spinmax,1)
                  else
                     spinmax = 1
                  endif
               else 
                  spinmax = 1
               endif

               do 230 spin = 1, spinmax
                  call soilbgc(spin,kpti,kptj)
230            continue

#ifndef SINGLE_POINT_MODEL
! determine the length of a precipitation event (between 4 and 24 hours),
! and time to start and end precipitation event. plen is in timesteps, while
! plens, startp, and endp are in seconds
               plenmin = 1 + int((4.0 * 3600. - 1.) / dtime)
               plenmax = max(int(24.0 * 3600. / dtime), plenmin)
        if(isimagro .eq. 0)then
           plen = min(plenmax,int(plenmin+ran2(seed,seed2,seed3,seed4) * (plenmax-plenmin+1)))
        else
           if ((iyear.lt.imetyear.or.iyear.gt.imetend.or.       &
              (iyear.eq.imetyear.and.jday.lt.dmetyear).or.       &
              (iyear.eq.imetend.and.jday.gt.dmetend)) ) then
              plen = min(plenmax,int(plenmin+ran2(seed,seed2,seed3,seed4) * (plenmax-plenmin+1)))
           else
              plen = 1
           endif
        endif ! check for crop existence
               plens   = dtime * plen
               startp  = dtime * min(niter-plen, int(ran2(seed,seed2,seed3,seed4)*(niter-plen+1)))
               endp    = startp + plens

#endif /* SINGLE_POINT_MODEL */

      if(isimagro .gt. 0) then
! calculate irrigation amount, timing, and duration
! if applicable to model run
!
! assume irrigation duration for a day is 12 hours long between 
! 6 am and 6 pm local time - this might be changed if managed irrigation can
! take place at night (e.g. crops)  
        ilens  = dtime * (12.0 * 3600. / dtime)   
        starti = dtime * (6.0  * 3600. / dtime)
        endi   = starti + ilens 
!
! only call irrigation if turned on and crops are planted
!
        if (irrigate .eq. 1) then 
           call irrigation(ilens,irrigate) 
        endif
      endif ! check for crop existence

!---------------------- HOUR HOUR HOUR -------------------------------------
! start of hourly loop
! if SINGLE_POINT_MODEL (0D), then this is half-hourly loop.
!---------------------- HOUR HOUR HOUR -------------------------------------
               do 240 istep = 1, niter

! calculate the time of day since midnight (in seconds)
                  if (dtime - int(dtime) .ne. 0) then
                     write(*,*) "ERROR: Time slice is not integer!"
                     stop 5
                  endif
                  timed = (istep - 1) * int(dtime)

#ifdef SINGLE_POINT_MODEL
! determine climatic conditions for the given timestep
                  call diurnal(timed, jday, test, dimforc)
! call the land surface model
                  call lsxmain(loopi, kpti, kptj, test, dimforc, isoilforc)
#else /* SINGLE_POINT_MODEL */
! determine climatic conditions for the given timestep

       if ((iyear.eq.imetyear.and.jday.ge.dmetyear).or.          &
            (iyear.gt.imetyear.and.iyear.lt.imetend ).or.        &
            (iyear.eq.imetend.and.jday.le.dmetend))then
                call methourly (iyear,jday,imonth,timed,seed, seed2,seed3,seed4)
                call diurnalmet(timed, jday, plens, startp, endp,  &
                     irrigate, ilens, starti, endi)
       else
                call diurnal(timed, jday, plens, startp, endp, starti,endi,seed, seed2,seed3,seed4)
       endif
! call the land surface model

! for now openmp is here only, because 80% of work is in lsxmain, 
! and diurnal cannot be parallelized without "fixing" weather generator
#ifdef _OPENMP
#ifndef SINGLE_POINT_MODEL
      !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(loopi,kpti,kptj,ithread,mpt)
      ithread = OMP_GET_THREAD_NUM()+1
      loopi = ithread
      kpti = 1 + (ithread-1)*(npoi/numthreads) 
      kptj = kpti + (npoi/numthreads) - 1  ! test when 1 numthreads
      if ( ithread .eq. numthreads ) then
         if ( kptj .ne. npoi ) then
            kptj = npoi
         end if
      end if 
      mpt = kptj - kpti + 1
      !print *,'ithread=',ithread,'thread=',OMP_GET_THREAD_NUM(),'kpti=',kpti,'kptj=',kptj
#endif
#endif

         if ( env_fastexec .lt. 2 ) then
                  call lsxmain(loopi, kpti, kptj)
         end if

#ifdef _OPENMP
#ifndef SINGLE_POINT_MODEL
      !$OMP END PARALLEL
#endif
#endif

#endif /* SINGLE_POINT_MODEL */

! call CTEM fire dynamics model
                   if ( isimveg.gt.0 .and. isimfire.eq.2 ) call firectem(kpti, kptj)

! accumulate some variables every timestep
                  call sumnow(kpti,kptj)
                  call sumday  (timed, loopi, kpti, kptj)
                  call summonth(timed, loopi, kpti, kptj)
                  call sumyear (timed, loopi, kpti, kptj)

! call to nitrogen stress routine for crops
               if(isimagro .gt. 0) then
                  call nitrostress(istep,iday)
!
                  call leaching(irestart,iyrrestart,istep,iholdsoiln)
               endif
! write out diagnostics
                  if (idiag .ne. 0) then
                     do 250 i = 1, idiag
                        if (iyear.ge.diagstart(i) .and. iyear.le.diagend(i) &
                            .and. ndiagpt(i).ge.1) then  
                           if (mod(istep,nfreq(i)).eq.0) then
                              call wdiag(i, istep)
                           end if
                        end if
250                  continue
                  end if

#ifndef SINGLE_POINT_MODEL
!this part is for calculating the crop 0D
    if(npoi .eq. 1 .and. isimagro .gt. 0)then
         call single_agro
    endif
#endif

#ifdef SINGLE_POINT_MODEL
! FIXME ET: this should not be inside the main loop, but in a function or subroutine...
! energy budget of the surface
         call single(linenum, test)
         if (test.eq.0.) goto 9999
#endif /* SINGLE_POINT_MODEL */

! ---------------------------------------------------------------
!               e n d   o f   t i m e   l o o p s
! ---------------------------------------------------------------

!---------------------- HOUR HOUR HOUR -------------------------------------
! end of the hourly loop
!---------------------- HOUR HOUR HOUR -------------------------------------
240            continue

!this part is for calculating the crop 0D

     if (isimagro .gt. 0) then
	   
         flx(10) = flx(10)/24. 
	     flx(14) = flx(14)/24. 
         flx(15) = flx(15)/24. 
 

!	     if((iyear.eq.2009.and.jday.gt.300).or.(iyear.eq.2010.and.jday.lt.140)) then
	     if(imetyear .ne. 9999) then

	              write(41,9201)iyear,jday,flx(1),flx(2),flx(3),flx(4),              &
                                flx(5),flx(6),flx(7),flx(8),flx(9),flx(11),flx(10),  &
                                plai(1,16)*grnfraccrop(1,16),plai(1,16),flx(12),     &
                                flx(13),flx(14),flx(15)
	endif  

 9201         format (2(i4,'  '),2(f6.2,'  '),3(f7.1,'  '),6(f6.2,'  '),6(f5.2,'  ')) 
     
        flx(1)= 0.0 
        flx(2)= 0.0
        flx(3)= 0.0
        flx(4)= 0.0
        flx(5)= 0.0
        flx(6)= 0.0
        flx(7)= 0.0
        flx(8)= 0.0
        flx(9)= 0.0
        flx(10)= 0.0     
        flx(11)= 0.0     
        flx(12)= 0.0     
        flx(13)= 0.0     
        flx(14)= 0.0     
        flx(15)= 0.0    
 
     endif! check for crop existence

#ifndef SINGLE_POINT_MODEL
! write out daily output
               if ( idailyout.eq.1 ) then
                  
                  if ( env_debug.gt.0  ) then
                     write (STDOUT,9101,advance='no') char(13), char(13), iyear, &
                          imonth, iday,              &
                          ' *writing daily output*'
#ifdef IFORT_NOFLUSH_BUG
                     close(STDOUT)
                     open(STDOUT)
#endif /* IFORT_NOFLUSH_BUG */
#ifdef PGF90_NOFLUSH_BUG
                     call flush(STDOUT)
#endif /* PFG90_NOFLUSH_BUG */
                  end if
!this part is for calculating the crop 0D
    if(npoi .gt. 1)then
                  call wdaily(jday,nday) 
    endif
               endif

! for long runs, write current day to file
         open(12,file=trim(outdir)//'/daysrun.dat',status='unknown')
         write(12,*) iyear,imonth,iday
         close(12)

#endif /* SINGLE_POINT_MODEL */
!------------------------- DAY DAY DAY -------------------------------------
! end of the daily loop
!------------------------- DAY DAY DAY -------------------------------------
220         continue
!
! write out monthly output (use 1st day of this month because of GrADS bug)
#ifndef SINGLE_POINT_MODEL
            if ( imonthout.eq.1 ) then

               if ( env_debug.gt.0 ) then
                  write (STDOUT,9101,advance='no') char(13), char(13), iyear, &
                       imonth, iday-1,            &
                       ' *writing monthly output*'
#ifdef IFORT_NOFLUSH_BUG
                  close(STDOUT)
                  open(STDOUT)
#endif /* IFORT_NOFLUSH_BUG */
#ifdef PGF90_NOFLUSH_BUG
                  call flush(STDOUT)
#endif /* PGF90_NOFLUSH_BUG */
               end if
!this part is not for calculating the crop 0D
    if(npoi .gt. 1)then
               call wmonthly(nday-ndaypm(imonth)+1)
    endif
            endif
#endif /* SINGLE_POINT_MODEL */
!---------------------- MONTH MONTH MONTH ----------------------------------
! end of the monthly loop
!---------------------- MONTH MONTH MONTH ----------------------------------
210      continue

#ifndef SINGLE_POINT_MODEL
! recalculate bioclimatic parameters if using anomalies
         if (iyear.ge.iyrdaily) call climanl2
#endif /* SINGLE_POINT_MODEL */

! perform vegetation dynamics
         if (isimveg.ne.0) call dynaveg(iwest,jnorth)

! calculate simple annual diagnostic variables for the globe
         call gdiag(snorth, ssouth, swest, seast, nrun, iyrlast, co2init, o2init, soilcspin)

! calculate simple annual diagnostic variables for each vegetation type
         call vdiag

#ifndef SINGLE_POINT_MODEL
! write out annual output
         if ( iyearout.eq.1 ) then

            if ( env_debug.gt.0 ) then
               write (STDOUT,9101,advance='no') char(13),char(13),iyear,imonth-1, &
                    iday-1, ' *writing yearly output*'
#ifdef IFORT_NOFLUSH_BUG
               close(STDOUT)
               open(STDOUT)
#endif /* IFORT_NOFLUSH_BUG */
#ifdef PGF90_NOFLUSH_BUG
               call flush(STDOUT)
#endif /* PGF90_NOFLUSH_BUG */
            end if
!this part is not for calculating the crop 0D
    if(npoi .gt. 1)then
! write yearly files
            call wyearly(nday)

! write HR dataset if requested
            if ( trim(hrmapfile) .ne. '' ) then
#ifdef HRMAP
               if ( env_fastexec .eq. 0 ) call whrvegmap(nday)
#else
               write(*,*) 'ERROR! hrmapfile defined in infile, but model not configured with --enable-highres_map option'
#endif
            endif

         endif
    endif

! write restart files
! Now just restart for dinamic vegetation (option 1 only)
         if(isimveg.ge.1) then

            if ( env_debug.gt.0 ) then
               write (STDOUT,9101,advance='no') char(13), char(13), iyear, imonth-1, &
                    iday-1, ' *writing restart data*'
#ifdef IFORT_NOFLUSH_BUG
               close(STDOUT)
               open(STDOUT)
#endif /* IFORT_NOFLUSH_BUG */
#ifdef PGF90_NOFLUSH_BUG
               call flush(STDOUT)
#endif /* PGF90_NOFLUSH_BUG */
            end if

! Now just restart for dinamic vegetation (option 1 only)
!this part is for calculating the crop 0D
    if(npoi .gt. 1)then
        if(isimveg.ge.1) call wrestart(nday, imonthlast)
    endif
        end if ! restart

#endif /* SINGLE_POINT_MODEL */

! update iyrlast value and file yearsrun.dat
         iyrlast = iyrlast + 1
!         open(12,file=trim(outdir)//'/yearsrun.dat',status='unknown')
!         write(filen,'(A16,I4,A4)')trim(outdir)//'/yearsrun.',iyrlast,'.dat'
         write(chyear(1:4),'(i4)') iyrlast
         filen = trim(outdir)//'/yearsrun.'//chyear//'.dat'
         open(12,file=filen,status='unknown')
         write(12,*) iyrlast
! uptade variables to turn possible restart of ran2 subroutine
         write(12,*) seed
         write(12,*) seed2
         write(12,*) seed3
         write(12,*) seed4
         write(12,*) aleafi
         write(12,*) astem
         write(12,*) astemi
         write(12,*) arepr
         write(12,*) iday
         write(12,*) cnleaf
         write(12,*) cnwood
         write(12,*) cnroot
         write(12,*) h20
         write(12,*) cndepth
         write(12,*) vf
         write(12,*) iniday
         write(12,*) fleaf
         write(12,*) fallrsgc
         close(12)
         write (STDOUT,9102) char(13), char(13), iyear, imonth-1, iday-1
!---------------------- YEAR YEAR YEAR -------------------------------------
! end of the yearly loop
!---------------------- YEAR YEAR YEAR -------------------------------------
200   continue

#ifdef SINGLE_POINT_MODEL
9999  continue
! If jumped here due to test = 0, then force calling the final g/vdiag.
      if (test.eq.0) then
         call gdiag
         call vdiag
      endif
#endif /* SINGLE_POINT_MODEL */

! end of the simulation
      write (*,*) ' '
      write (*,*) '*** end of run ***'

! timing report - only works with gfortran
#ifdef __GFORTRAN__
      call ETIME(tarray, tresult)  
      telapsed = time() - telapsed
      write(*,*) ' '
      write(*,1001) 'elapsed time:', INT(telapsed), REAL(telapsed)/60.0
      write(*,1001) 'run time:    ', INT(tresult), tresult/60.0
      write(*,1001) 'user time:   ', INT(tarray(1)), tarray(1)/60.0
      write(*,1001) 'system time: ', INT(tarray(2)), tarray(2)/60.0
      write(*,*) ' '
1001           format (a,i10,' s (',f10.1,' m)')
1002           format (a,i4,a,i4,a,i4,'s (',f3.1,' m)')
#endif

      stop

end program main !main_offline
