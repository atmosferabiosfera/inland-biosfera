#include "inland_config.h"
#include "inland_compar.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
subroutine wrestart (nday, imonthp)
! ---------------------------------------------------------------------
! instantaneous output for restarts
!---------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_control, only: iyear, iyear0, outdir, env_floatout, env_debug
      use inland_comveg
      use inland_comsoi
      use inland_comsno
      use inland_comwork
      use inland_comsum
      use inland_comatm, only: xstore
      use inland_combcs, only: xinwet
      use inland_comcrop
      use inland_comnitr

      implicit none
!-----------------------------------------------------------------------
! this subroutine writes the restart values of:
!
!  fi      = fractional snow cover
!  tsno    = temperature of snow
!  hsno    = snow depth
!  tsoi    = soil temperature
!  wisoi   = soil ice content
!  wsoi    = soil moisture content
!  cbiol   = carbon in leaf biomass pool
!  cbiow   = carbon in woody biomass pool
!  cbior   = carbon in fine root biomass pool
!  sapfrac = sapwood fraction
!  clitlm  = leaf metabolic litter
!  clitls  = leaf structural litter
!  clitll  = leaf lignin litter
!  clitrm  = root metabolic litter
!  clitrs  = root structural litter
!  clitrl  = root lignin litter
!  clitwm  = woody metabolic litter
!  clitws  = woody structural litter
!  clitwl  = woody lignin litter
!  falll   = annual leaf litterfall
!  fallr   = annual fine root turnover 
!  fallw   = annual wood litterfall
!  fu      = fraction of overall area covered by upper canopy
!  fl      = fraction of snow-free area covered by lower  canopy
!  vegtype0= vegetation type
!  totcmic = total microbial carbon
!  csoislop= slow soil carbon, protected humus
!  csoislon= slow soil carbon, nonprotected humus
!  csoipas = passive soil carbon
!  gdd0    = growing degree days 0
!  gdd5    = growing degree days 5
!  tc      = coldest monthly temperature
!  tw      = warmest monthly temperature
!  gdd0this= growing degree days 0 of the current year
!  gdd5this= growing degree days 5 of the current year
!  tcthis  = coldest monthly temperature of the current year
!  twthis  = warmest monthly temperature of the current year
!  wipud   = ice content of puddles per soil area
!  wpud    = liquid content of puddles per soil area
!  wliql   = intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
!  wliqs   = intercepted liquid h2o on upper canopy stem area (kg m-2)
!  wliqu   = intercepted liquid h2o on upper canopy leaf area (kg m-2)
!  wsnol   = intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
!  wsnos   = intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
!  wsnou   = intercepted frozen h2o (snow) on upper canopy leaf area
!     (kg m-2)
!  tu      = temperature of upper canopy leaves (K)
!  ts      = temperature of upper canopy stems (K)
!  tl      = temperature of lower canopy leaves & stems(K)
!  t12     = air temperature at z12 (K)
!  t34     = air temperature at z34 (K)
!  q12     = specific humidity of air at z12
!  q34     = specific humidity of air at z34
!
!
!  fractional cover of upper and lower canopies
!  
!  gdd8clim    = growing degree days > 8C for corn crops between April 1 - Sept 30
!  gdd10clim   = growing degree days > 10C for soybean crops between April 1 - Sept 30
!  gdd12clim   = growing degree days > 12C for sugarcane crops between April 1 - Sept 30
!  gdd0clim    = growing degree days > 0C
!  currgdd8    = annual total growing degree days for current year
!  currgdd10   = annual total growing degree days for current year
!  currgdd12   = annual total growing degree days for current year
!  currgdd0    = annual total growing degree days for current year
!  nfixnat     = annual total nitrogen fixation by natural vegetation (kg n m-2 y-1)
!  ftot   = annual total inorganic nitrogen leached from soil profile (kg n  m-2 y-1)
!  no3leach    = annual total nitrate-n leached from profile (kg n m-2 y-1)
!  nbalance    = annual soil nitrogen balance calculation (kg n ha-1)
!  natvegnup   = annual total inorganic nitrogen uptake by natural vegetation (kg n m-2 y-1)
!  decomps     = soil organic matter decomposition factor     (dimensionless)
!  decompl     = litter decomposition factor                  (dimensionless)



!  Physiology
!  ciub    = intercellular co2 concentration - broadleaf (mol_co2/mol_air)
!  ciuc    = intercellular co2 concentration - conifer   (mol_co2/mol_air)
!  cils    = intercellular co2 concentration - shrubs    (mol_co2/mol_air)
!  cil3    = intercellular co2 concentration - c3 plants (mol_co2/mol_air)
!  cil4    = intercellular co2 concentration - c4 plants (mol_co2/mol_air)
!  csub    = leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
!  csuc    = leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
!  csls    = leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
!  csl3    = leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
!  csl4    = leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
!  gsub    = upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
!  gsuc    = upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
!  gsls    = lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
!  gsl3    = lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
!  gsl4    = lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
!
!  phenology :
!  agddu   = annual accumulated growing degree days for bud burst, upper canopy
!  agddl   = annual accumulated growing degree days for bud burst, lower canopy
!  tempu   = cold-phenology trigger for trees
!  templ   = cold-phenology trigger for grasses/shrubs
!  td      = daily average air temperature
!  a10td    = 10-day avg daily temp
!  a10ancub = 10-day average canopy photosynthesis rate - broadleaf
!  a10ancuc = 10-day average canopy photosynthesis rate - conifers
!  a10ancls = 10-day average canopy photosynthesis rate - shrubs
!  a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
!  a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
!  a10scalparamu = 10-day average canopy scaling parameter - upper canopy
!  a10scalparaml = 10-day average canopy scaling parameter - lower canopy
!  a10daylightu = 10-day average daylight - upper canopy
!  a10daylightl = 10-day average daylight - lower canopy
!  dropu   = drought-phenology trigger for trees
!  dropls  = drought-phenology trigger for shrubs
!  dropl4  = drought-phenology trigger for c4 grasses
!  dropl3  = drought-phenology trigger for c3 grasses
!  lail = lai(i,1) : lower canopy lai
!  laiu = lai(i,2) : upper canopy lai
! (NOTE: a10ancuc is not used at this point, so its wrestart entry 
! is commented out)
!
! Summation variables used in dynaveg, gdiag and vdiag
!
!  aynpp = annual total npp per pft (kg-c/m**2/yr)
!  aygpp = annual total gpp per pft (kg-c/m**2/yr)
!  ayalit = yearly av. aboveground litter (kg-c/m**2)
!  ayblit = yearly av. belowground litter (kg-c/m**2)
!  aycsoi = yearly av. soil carbon content (kg-c/m**2)
!  ayco2root = yearly av. root respiration (kg-C/m**2/yr)
!  ayco2mic = yearly av. microbial respiration (kg-C/m**2/yr)
!  yrleach = yearly av. C leaching rate (kg_C m-2/yr)

!  ayanlit = yearly av. aboveground litter nitrogen (kg-N/m**2)
!  aybnlit = yearly av. belowground litter nitrogen (kg-N/m**2)
!  aynsoi = yearly av. total soil nitrogen (kg-N/m**2)
!  ynleach = yearly av. N leaching rate (kg_C m-2/yr)

!  ayprcp = yearly av. precipitation (mm/year)
!  ayaet  = yearly av. aet (mm/yr)
!  aytrans = yearly av. transpiration (mm/yr)
!  aysrunoff = yearly av. surface runoff (mm/yr)
!  aydrainage = yearly av. drainage (mm/yr)
!  dwtot = yearly h20 change (mm/yr)
!  wtot = yearly h20 content (mm/yr)
!
!  aysolar = annual average incident solar radiation (w/m**2)
!  ayreflect = annual average reflected radiation (W/m**2) 
!  ayirup = annual average upward ir radiation (w/m**2)
!  ayirdown = annual average downward ir radiation (w/m**2)
!  aysens = annual average sensible heat flux (w/m**2)
!  aylatent = annual average latent heat flux (w/m**2)
!
!  aywsoi = yearly av. soil moisture content
!  aywisoi = yearly av. soil ice content
!  aytsoi = yearly av.  temperature
!  ayvwc = annual average 1m volumetric water content (fraction)
!  ayawc = annual average 1m plant-available water content (fraction)
!
!  firefac = annual average fire factor
!  ayrootbio = annual average live root biomass (kg-C/m**2/yr)
!  aynmintot = annual total nitrogen mineralization (kg-N/m**2/yr)
!  aycmic = annual av. microbial biomass (kg-c/m**2)

!Crops
!
!  aplantn = available inorganic nitrogen in soil profile for plant growth (kg_n m-2 y-1)
!  smsoil = current timestep solute in soil (kg_solute/m2)
!  smsoln = current timestep solute in solution (kg_solute/m2)
!  cbios = carbon in stem biomass pool (kg_C m-2)
!  aybprod = annual total biomass production for crops
!  ayabprod = annual total aboveground biomass production for crop
!  ayrprod = annual total root carbon accumulation for crops (kg-c/m**2/yr)
!  aylprod = annual total leaf carbon accumulation for crops (kg-c/m**2/yr)
!  harvdate = day of year that crop pft was harvested
!  cropplant = index keeping track of whether that crop has been planted during current year
!  biomass = biomass
!  hui = heat unit index
!  idop = day of year that crop was planted
!  fnleaf = current fraction of nitrogen in leaf dry matter
!  fnstem = current fraction of nitrogen in stem dry matter
!  stressn = current fraction of nitrogen in stem dry matter
!  thrlai = lai threshold for crops when senescence begins
!  grnfraccrop = green fraction of aboveground vegetation for crops
!  dtnleach = daily total inorganic nitrogen solute (ammonium and nitrate) leached
!  croplive = 0 crops have been planted and living : 1 crops not living
!  gddplant = accumulated growing degree days past planting date for crop = j
!  gddtsoi = accumulated growing degree days past planting date for crop = j
!  idpp = number of days past planting
!  pstart = day to start plant, relative to the planting calendar
!  cntops = cn ratio of plant residue
!  pdate = planting date (real value)
!  hdate = harvest date (real value)
!  harvidx = end of year harvest index for crop
!  croplaimx = maximum attained lai by crop during growing season
!  cropyld = crop yield in bu/ac
!  dmyield = yield dry matter in Mg/ha
!  dmleaf = leaf dry matter in Mg/ha
!  dmstem = stem dry matter in Mg/ha
!  dmroot = fine root dry matter (Mg/ha)
!  dmresidue = aboveground leaf and stem residue dry matter in Mg/ha
!  residuen = nitrogen contained in aboveground crop residue
!  nconcl = end of season N concentration in leaf (fraction)
!  nconcs = end of season N concentration in stem (fraction)
!  nconcr = end of season N concentration in root (fraction)
!  nconcg = end of season N concentration in grain (fraction)
!  cropn = nitrogen removed by crop in kg/ha/yr
!  cropfixn = nitrogen fixation by crop in kg/ha/yr
!  cnrootvec = cn ratio of plant roots (renamed to avoid conflict with inlands cnroot)
!  crmclim = crop relative maturity rating (CRM) from Pioneer regression relationships
!  crmact = best crop relative maturity rating (CRM) for that year based on total GDD accumulated
!  crmplant = crop relative maturity rating (CRM) based on GDD accumulated since planting of that year
!  gddmaturity = accumulated growing degrees needed for plant to reach both vegetative and physiological maturity
!  cropout(:,4,35) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,5,35) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,6,35) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,13,50) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,14,50) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,15,50) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,16,50) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,13,52) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,14,52) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,15,52) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,16,52) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,13,51) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,14,51) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,15,51) = crop cycle number in a cycle (1s planting - last ratoon)
!  cropout(:,16,51) = crop cycle number in a cycle (1s planting - last ratoon)
!  adnpp = daily total npp for each plant type (kg-C/m**2/day)
!  tnpp = total NPP per time step for each pft (kg_C m-2 / timestep)
!  plaimx = annual maximum leaf area index of each plant functional type
!  plai = initial total LAI for each vegtype (used in iniveg)
!  totnuptake = annual total nitrogen uptake (kg_n m-2 y-1)
!  totnfix = annual total nitrogen fixation through symbiosis (kg_n m-2 y-1)
!  cbiog = carbon in grain biomass pool (kg_C m-2)
!  gdd8 = growing degree days > 8C for corn crops
!  gdd10 = growing degree days > 10C for soybean crops
!  gdd12 = growing degree days > 12C for sugarcane crops
!  gdd0c = growing degree days > 0C for wheat crops
!  gdd0clim = climate gdd0
!  ftot = total leaching of inorganic  nitrogen
!  gddpl15 = accumulated growing degree days for 1 and half year for sugarcane
!  gddsgcp1 = factor to get the hybrid planted for sugarcane plating, 1 and half year
!  gddsgcp2 = factor to get the hybrid planted for sugarcane plating, 1 and half year
!  htmx1 = maximum height attained by a crop during year
!  htmx2 = maximum height attained by a crop during year
!  ik = counter for sugarcane crop one and half year GDD
!  cropy = sugarcane crop year since planting
!  cdays = days since planting in weather.f
!  ncyears = number of the ith crop year
!  ctot = total inorganic nitrogen in soil profile (kg n  m-2)
!  drntot = total inorganic nitrogen in soil profile (kg n  m-2)
!  cic3 = intercellular co2 concentration - c3 crops (mol_co2/mol_air)
!  cic4 = intercellular co2 concentration - c4 crops (mol_co2/mol_air)
!  csc3 = leaf boundary layer co2 concentration - c3 crops (mol_co2/mol_air)
!  csc4 = leaf boundary layer co2 concentration - c4 crops (mol_co2/mol_air)
!  ua = wind speed (m s-1)
!  gddemerg = accumulated growing degree days at leaf emergence

! input variables
      integer nday,    & ! number of days run since iyear0
              imonthp    ! this calendar month

! local variables
      integer nyears,  & ! years of run iyear0
              mstep, lf, lenchrib, n, idies, istat, k

      integer istart(5), &
              icount(5)   ! for writing restart vars

      character*21 tunits
      character*20 fdir         ! used to construct odd/even file names
      character*8 tmpdate      ! temp. variable to hold date
      character*8 cdate        ! date to use in history attribute in files
      character*80 dimnames(5)  ! names of dimensions for restart vars
      integer ndims        ! number of dims (including tile, not including pft)
      character*1024 filen      ! file name
      character*4  chyear
      character*2  chmonth

      real*8 slayers(nsnolay),  & ! index for snow layers
             depthsoi(nsoilay), & ! soil layer depths
             pindex(npft),      & ! index for pfts
             tindex(mlpt),      & ! index used for tiles
             xindex(3),         & ! index for xstore
             ftime                ! floating point time value
      
      integer tmp_floatout

! extra dim variables passed to inifile_dims
      integer, parameter :: numXdim = 5 ! snow, soil, pft, tile, xstore
      integer, parameter :: numXvals = 100 ! pick a big number to make sure we have enough space
      integer nXth(numXdim)
      character*80, dimension(numXDim) :: nameXth, longXth, unitsXth, axisXth
      real*8 valsXth(numXVals,numXdim)

      data istart / 1,1,1,1,1 /, &
           icount / nlon,nlat,1,1,1 /
! ---------------------------------------------------------------------

      ! disable float output for restart
      tmp_floatout = env_floatout
      env_floatout = 0

      icount(1) = nlonsub
      icount(2) = nlatsub

      mstep = 1
      chyear = '0000'
      if (iyear .gt. 1000) then
         write(chyear(1:4),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(chyear(4:4),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(chyear(3:4),'(i2)') iyear
      else
         write(chyear(2:4),'(i3)') iyear
      endif

      chmonth = '00'
      if (imonthp .gt. 9) then
         write(chmonth(1:2),'(i2)') imonthp
      else
         write(chmonth(2:2),'(i1)') imonthp
      endif

      fdir = 'restart/'
      lf = lenchrib(fdir)

! how many dims?
      if ( mlpt .gt. 1 ) then
         ndims = 4
      else
         ndims = 3
      end if

! initialize snow layer indicies, pft names, etc (done as in wyearly)
      if (mstep .eq. 1) then
         call date_and_time(tmpdate)
         cdate(1:2) = tmpdate(5:6)
         cdate(3:3) = '/'
         cdate(4:5) = tmpdate(7:8)
         cdate(6:6) = '/'
         cdate(7:8) = tmpdate(3:4)
         ftime = nday

! time units is days since Dec 31 of the year before iyear0
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1

         do n = 1, nsnolay
            slayers(n) = float(n)
         end do

         depthsoi(1) = hsoi(1)
         do n = 2, nsoilay
            depthsoi(n) = depthsoi(n-1)+hsoi(n)
         end do

! define other dimension indexes
         do n = 1, npft
            pindex(n) = n
         end do

         do n = 1, mlpt
            tindex(n) = n
         end do

         do n = 1, 3
            xindex(n) = float(n)
         end do

! define dimnames
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
! dimnames(3+) is set for each variable seperately
         dimnames(3:) = 'time'

! fill variables passed to inifile_dims

         ! pfts
         n = 1
         nameXth(n)  = 'pft'
         longXth(n)  = 'plant functional type'
         unitsXth(n) = 'none'
         axisXth(n)  = 'E'
         nXth(n)     = npft
         valsXth(1:npft,n) = pindex

         ! tiles
         n = 2
         nameXth(n)  = 'tile'
         longXth(n)  = 'tile'
         unitsXth(n) = 'none'
         axisXth(n)  = 'Z'
         nXth(n)     = mlpt
         valsXth(1:mlpt,n) = tindex

         ! xstore
         n = 3
         nameXth(n)  = 'xstorelevel'
         longXth(n)  = 'level'
         unitsXth(n) = 'none'
         axisXth(n)  = 'Z'
         nXth(n)     = 3
         ! up/down not used...
         valsXth(1:3,n) = xindex

         ! snow layers
         n = 4
         nameXth(n)  = 'snowlayer'
         longXth(n)  = 'snow layers top to bottom'
         unitsXth(n) = 'none'
         axisXth(n)  = 'Z'
         nXth(n)     = nsnolay
         ! up/down not used...
         valsXth(1:nsnolay,n) = slayers

         ! soil layers
         n = 5
         nameXth(n)  = 'soillayer'
         longXth(n)  = 'depth of soil layer bottom'
         unitsXth(n) = 'meter'
         axisXth(n)  = 'Z'
         nXth(n)     = nsoilay
         ! up/down not used...
         valsXth(1:nsoilay,n) = depthsoi

      endif

! by Etienne on 06/01/2013
! one restart file with all vars and dims!
      if (myid .eq. 0) then
         filen = trim(outdir)//'/inland-restart-'//chyear//'_'//chmonth//'.nc'

         if ( env_debug.gt.1 ) print *,'creating restart file '//trim(filen)

         if (mstep .eq. 1) then

!-----------------------------------------------------------------------
! inifile
            call inifile_dims(idies,filen,'restart file',   &
                              'inland wrestart',cdate,nlonsub,lonscale,nlatsub,latscale, &
                              nameXth,longXth,unitsXth,axisXth,nXth,valsXth,numXdim,numXVals, &
                              tunits,'gregorian',istat) 

!-----------------------------------------------------------------------
! inivar

! ===== xstore
            dimnames(3) = 'xstorelevel'
            !dimnames(4) = 'tile' ! for now xstore is from 1:npoi1, so no tile
            dimnames(4) = 'time'

            call inivar(idies,'xstore','xstore value', &
                        'none',4,dimnames,istat)

! ===== snow layers
            dimnames(3) = 'snowlayer'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'

            call inivar(idies,'tsno','instantaneous snow cover temperature', &
                        'degK',ndims+1,dimnames,istat)
            call inivar(idies,'hsno','instantaneous snow layer thickness', &
                        'meters',ndims+1,dimnames,istat)

! ===== soil layers
            dimnames(3) = 'soillayer'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'

            call inivar(idies,'tsoi','instantaneous soil temperature','degK', &
                        ndims+1,dimnames,istat)
            call inivar(idies,'wisoi','instantaneous fraction of soil pore '// &
                        'space containing ice','fraction',ndims+1,dimnames,  &
                        istat)
            call inivar(idies,'wsoi','instantaneous fraction of soil pore '// &
                        'space containing water','fraction',ndims+1,dimnames, &
                        istat)
!Crops
         if(isimagro .gt. 0)then
            call inivar(idies,'aplantn','available inorganic nitrogen in soil profile for plant growth (kg_n m-2 y-1)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'smsoil','current timestep solute in soil (kg_solute/m2)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'smsoln','current timestep solute in solution (kg_solute/m2)', &
                        '(days/month)',ndims+1,dimnames,istat)
         endif


! ===== Carbon pools and fluxes per pft
            dimnames(3) = 'pft'
            dimnames(4:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(4) = 'tile'

            call inivar(idies,'cbiol','instantaneous carbon in leaf biomass '//&
                        'pool','kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'cbiow','instantaneous carbon in wood biomass '//&
                        'pool','kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'cbior','instantaneous carbon in root biomass '//&
                        'pool','kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'cbios','carbon in stem biomass pool (kg_C m-2) '//&
                        'pool','kg/m^2',ndims+1,dimnames,istat)
            call inivar(idies,'aynpp','current value of annual sum of npp', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'aygpp','current value of annual sum of gpp', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
! ====Crops
         if(isimagro .gt. 0)then
            call inivar(idies,'aybprod','annual total biomass production for crops', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'ayabprod','annual total aboveground biomass production for crop', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'ayrprod','annual total root carbon accumulation for crops (kg-c/m**2/yr)', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'aylprod','annual total leaf carbon accumulation for crops (kg-c/m**2/yr)', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'harvdate','day of year that crop pft was harvested', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'cropplant','index keeping track of whether that crop has been planted during current year', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'biomass','biomass', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'hui','heat unit index', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'idop','day of year that crop was planted', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'fnleaf','current fraction of nitrogen in leaf dry matter', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'fnstem','current fraction of nitrogen in stem dry matter', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'stressn','current fraction of nitrogen in stem dry matter', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'thrlai','lai threshold for crops when senescence begins', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)

            call inivar(idies,'grnfraccrop','green fraction of aboveground vegetation for crops', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)
            call inivar(idies,'dtnleach','daily total inorganic nitrogen solute (ammonium and nitrate) leached', &
                        'kg-C/m^2/yr',ndims+1,dimnames,istat)

            call inivar(idies,'croplive','0 crops have been planted and living : 1 crops not living', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'gddplant','accumulated growing degree days past planting date for crop = j', &
                        '(days/month)',ndims+1,dimnames,istat)

            call inivar(idies,'gddtsoi','accumulated growing degree days past planting date for crop = j', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'idpp','number of days past planting', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'pstart','day to start plant, relative to the planting calendar', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'cntops','cn ratio of plant residue', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'pdate','planting date (real value)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'hdate','harvest date (real value)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'harvidx','end of year harvest index for crop', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'croplaimx','maximum attained lai by crop during growing season', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'cropyld','crop yield in bu/ac', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'dmyield','yield dry matter in Mg/ha', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'dmleaf','leaf dry matter in Mg/ha', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'dmstem','stem dry matter in Mg/ha', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'dmroot','fine root dry matter (Mg/ha)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'dmresidue','aboveground leaf and stem residue dry matter in Mg/ha', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'residuen','nitrogen contained in aboveground crop residue', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'nconcl','end of season N concentration in leaf (fraction)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'nconcs','end of season N concentration in stem (fraction)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'nconcr','end of season N concentration in root (fraction)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'nconcg','end of season N concentration in grain (fraction)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'cropn','nitrogen removed by crop in kg/ha/yr', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'cropfixn','nitrogen fixation by crop in kg/ha/yr', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'cnrootvec','cn ratio of plant roots (renamed to avoid conflict with inlands cnroot)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'crmclim','crop relative maturity rating (CRM) from Pioneer regression relationships', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'crmact','best crop relative maturity rating (CRM) for that year based on total GDD accumulated', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'crmplant','crop relative maturity rating (CRM) based on GDD accumulated since planting of that year', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'gddmaturity','accumulated growing degrees needed for plant to reach both vegetative and physiological maturity', &
                        '(days/month)',ndims+1,dimnames,istat)

            call inivar(idies,'cropout(:,4,35)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,5,35)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,6,35)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,13:,50)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,14:,50)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,15:,50)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,16:,50)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)

            call inivar(idies,'cropout(:,13:,52)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,14:,52)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,15:,52)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,16:,52)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)

            call inivar(idies,'cropout(:,13:,51)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,14:,51)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,15:,51)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)
            call inivar(idies,'cropout(:,16:,51)','crop cycle number in a cycle (1s planting - last ratoon)', &
                        'day of year',ndims+1,dimnames,istat)


            call inivar(idies,'adnpp','daily total npp for each plant type (kg-C/m**2/day)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'tnpp','total NPP per time step for each pft (kg_C m-2 / timestep)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'plaimx','annual maximum leaf area index of each plant functional type', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'plai','initial total LAI for each vegtype (used in iniveg)', &
                        '(days/month)',ndims+1,dimnames,istat)

            call inivar(idies,'totnuptake','annual total nitrogen uptake (kg_n m-2 y-1)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'totnfix','annual total nitrogen fixation through symbiosis (kg_n m-2 y-1)', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'cbiog','carbon in grain biomass pool (kg_C m-2)', &
                        '(days/month)',ndims+1,dimnames,istat)

            call inivar(idies,'aroot','carbon allocation fraction to fine roots', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'aleaf','carbon allocation fraction to leaves', &
                        '(days/month)',ndims+1,dimnames,istat)
            call inivar(idies,'awood','carbon allocation fraction to wood', &
                        '(days/month)',ndims+1,dimnames,istat)
        endif
            
! All other 2-d variables
            dimnames(3:) = 'time'
            if ( mlpt .gt. 1 ) dimnames(3) = 'tile'

            call inivar(idies,'anytime','intrayear timestep index','dimless',  &
                        ndims,dimnames,istat)
            call inivar(idies,'fi','instantaneous fractional snow cover', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'sapfrac','instantaneous sapwood fraction', &
                        'fraction',ndims,dimnames,istat)
            call inivar(idies,'clitlm','instantaneous leaf metabolic litter'// &
                        'carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitls','instantaneous leaf structural litter'//&
                        'carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitll','instantaneous leaf lignin litter '// &
                        'carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitrm','instantaneous root metabolic litter '//&
                        'carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitrs','instantaneous root structural litter'//&
                        ' carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitrl','instantaneous root lignin litter '// &
                        'carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitwm','instantaneous woody metabolic litter'//&
                        ' carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitws','instantaneous woody structural '// &
                        'litter carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'clitwl','instantaneous woody lignin litter '// &
                        'carbon','kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'falll','annual leaf litterfall carbon','kg/m^2',&
                        ndims,dimnames,istat)
            call inivar(idies,'fallr','annual fine root turnover carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'fallw','annual wood turnover carbon','kg/m^2', &
                        ndims,dimnames,istat)
            call inivar(idies,'vegtype0','vegetation type','none',ndims,dimnames, &
                        istat)
            call inivar(idies,'fu','upper canopy fractional cover','none',ndims, &
                        dimnames,istat)
            call inivar(idies,'fl', 'upper canopy fractional cover (snow-free area)', &
                        'none',ndims,dimnames,istat)

         if(isimagro .gt. 0)then
            call inivar(idies,'gdd8','growing degree days > 8C for corn crops', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd10','growing degree days > 10C for soybean crops', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd12','growing degree days > 12C for sugarcane crops', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd0c','growing degree days > 0C for wheat crops', &
                        'degrees C',ndims,dimnames,istat)
            call inivar(idies,'gdd0clim','climate gdd0', &
                        'degrees C',ndims,dimnames,istat)

            call inivar(idies,'ftot','annual total inorganic nitrogen leached from soil profile (kg n  m-2 y-1)', &
                       'kg/hectare',ndims,dimnames,istat)
         endif
            call inivar(idies,'totcmic','instantaneous total microbial carbon',&
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'csoislop','instantaneous slow soil carbon protected humus', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'csoislon', 'instantaneous slow soil carbon nonprotected humus', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'csoipas','instantaneous passive soil carbon', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'decomps','soil decomposition factor','none',ndims, &
                        dimnames,istat)
            call inivar(idies,'decompl','litter decomposition factor','none', &
                        ndims,dimnames,istat)

! TODO: theoretically fixed vegetation runs does not require gdd0/5 gddthis0/5
!      AND pfts information to be written. Then we should disable writing this
!      on static vegetation runs!
            call inivar(idies,'gdd0','instantaneous growing degree days above 0 deg_C', &
                        'days degC',ndims,dimnames,istat)
! TODO: theoretically fixed vegetation runs does not require gdd0/5 gddthis0/5
!      AND pfts information to be written. Then we should disable writing this
!      on static vegetation runs!
            call inivar(idies,'gdd5','instantaneous growing degree days above 5 deg_C', &
                        'days degC',ndims,dimnames,istat)
            call inivar(idies,'tc','instantaneous coldest monthly temperature',&
                        'degC',ndims,dimnames,istat)
            call inivar(idies,'tw','instantaneous warmest monthly temperature',&
                        'degC',ndims,dimnames,istat)
! TODO: theoretically fixed vegetation runs does not require gdd0/5 gddthis0/5
!      AND pfts information to be written. Then we should disable writing this
!      on static vegetation runs!
            call inivar(idies,'gdd0this','instantaneous growing degree days above 0 deg_C for current year', &
! TODO: theoretically fixed vegetation runs does not require gdd0/5 gddthis0/5
!      AND pfts information to be written. Then we should disable writing this
!      on static vegetation runs!
                        'days degC',ndims,dimnames,istat)
            call inivar(idies,'gdd5this','instantaneous growing degree days above 5 deg_C for current year', &
                        'days degC',ndims,dimnames,istat)
            call inivar(idies,'tcthis','instantaneous coldest monthly temperature for current year', &
                        'degC',ndims,dimnames,istat)
            call inivar(idies,'twthis','instantaneous warmest monthly temperature for current year', &
                        'degC',ndims,dimnames,istat)

            call inivar(idies,'wipud','instantaneous ice content of puddles', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'wpud','instantaneous liquid water content of puddles', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'wliql', &
                        'instantaneous intercepted water on lower canopy', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'wliqs', &
                        'instantaneous intercepted water on upper canopy stems', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'wliqu', &
                        'instantaneous intercepted water on upper canopy leaves', &
                        'kg/m^2',ndims,dimnames,istat)

            call inivar(idies,'wsnol', &
                        'instantaneous intercepted snow on lower canopy', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'wsnos', &
                        'instantaneous intercepted snow on upper canopy stems', &
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'wsnou', &
                       'instantaneous intercepted snow on upper canopy leaves',&
                        'kg/m^2',ndims,dimnames,istat)
            call inivar(idies,'tu', &
                        'instantaneous temperature of upper canopy leaves', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'ts', &
                        'instantaneous temperature of upper canopy stems', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'tl', &
                        'instantaneous temperature of lower canopy leaves', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'t12', &
                        'instantaneous temperature at level z12', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'t34', &
                        'instantaneous temperature at level z34', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'tlsub', &
                      'instantaneous temperature of snow burried lower canopy',&
                        'K',ndims,dimnames,istat)
            call inivar(idies,'tg', &
                        'instantaneous soil skin temperature', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'ti', &
                        'instantaneous snow skin temperature', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'q12', &
                        'instantaneous humidity at level z12', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'q34', &
                        'instantaneous humidity at level z34', &
                        'K',ndims,dimnames,istat)
            call inivar(idies,'ciub', &
                        'intercellular co2 concentration - broadleaf', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'ciuc', &
                        'intercellular co2 concentration - conifer', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'cils', &
                        'intercellular co2 concentration - shrub', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'cil3', &
                        'intercellular co2 concentration - c3 grasses', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'cil4', &
                        'intercellular co2 concentration - c4 grasses', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'csub', &
                        'leaf boundary layer co2 concentration - broadleaf', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'csuc', &
                        'leaf boundary layer co2 concentration - conifer', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'csls', &
                        'leaf boundary layer  co2 concentration - shrub', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'csl3', &
                        'leaf boundary layer  co2 concentration - c3 grasses', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'csl4', &
                        'leaf boundary layer co2 concentration - c4 grasse', &
                        'mol_co2/mol_air',ndims,dimnames,istat)
            call inivar(idies,'gsub', &
                        'upper canopy stomatal conductance - broadleaf', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'gsuc', &
                        'upper canopy stomatal conductance - conifer', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'gsls', &
                        'upper canopy stomatal conductance - shrub', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'gsl3', &
                        'upper canopy stomatal conductance - c3 grasses', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'gsl4', &
                        'upper canopy stomatal conductance - c4 grasses', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'agddu', &
                        'instantaneous growing degree days', &
                        'days degC',ndims,dimnames,istat)
            call inivar(idies,'agddl', 'instantaneous ','days degC',ndims,dimnames,&
                        istat)
            call inivar(idies,'tempu','cold phenology trigger for trees', &
                        'dimensionless',ndims,dimnames,istat)
            call inivar(idies,'templ', &
                        'cold phenology trigger for grasses/shrubs', &
                        'dimensionless',ndims,dimnames,istat)
            call inivar(idies,'td','average daily air T','K',ndims,dimnames, &
                        istat)
            call inivar(idies,'a10td','10-day average daily air T','K',ndims, &
                        dimnames,istat)
            call inivar(idies,'a10ancub', &
                        '10-day average canopy photosynth. rate, broadleaf', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'a10ancls', &
                        '10-day average canopy photosynth. rate, shrubs', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'a10ancl4', &
                        '10-day average canopy photosynth. rate, c4 grasses', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'a10ancl3', &
                        '10-day average canopy photosynth. rate, c3 grasses', &
                        'mol_co2 m-2 s-1',ndims,dimnames,istat)
            call inivar(idies,'a10scalparamu','10-day average canopy scaling parameter, upper canopy', &
                        'none',ndims,dimnames,istat)
            call inivar(idies,'a10scalparaml','10-day average canopy scaling parameter, lower canopy', &
                        'none',ndims,dimnames,istat)
            call inivar(idies,'a10daylightu','10-day average daylight, upper canopy', &
                        'W m-2',ndims,dimnames,istat)
            call inivar(idies,'a10daylightl','10-day average daylight, lower canopy', &
                        'W m-2',ndims,dimnames,istat)
            call inivar(idies,'dropu','drought-pheno. trigger for trees', &
                        'dimensionless',ndims,dimnames,istat)
            call inivar(idies,'dropls','drought-pheno. trigger for shrubs', &
                        'dimensionless',ndims,dimnames,istat)
            call inivar(idies,'dropl4','drought-pheno. trigger for c4 grasses', &
                        'dimensionless',ndims,dimnames,istat)
            call inivar(idies,'dropl3','drought-pheno. trigger for c3 grasses', &
                        'dimensionless',ndims,dimnames,istat)
            call inivar(idies,'laiu','LAI of upper canopy','dimensionless',ndims, &
                        dimnames,istat)
            call inivar(idies,'lail','LAI of lower canopy','dimensionless',ndims, &
                        dimnames,istat)


            call inivar(idies,'ayalit','yearly av. aboveground litter', &
                        '(kg-c/m**2)',ndims,dimnames,istat)
            call inivar(idies,'ayblit','yearly av. belowground litter', &
                        '(kg-c/m**2)',ndims,dimnames,istat)
            call inivar(idies,'aycsoi','yearly av. soil carbon content', &
                        '(kg-c/m**2)',ndims,dimnames,istat)
            call inivar(idies,'ayco2root','yearly av. root respiration', &
                        '(kg-C/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'ayco2mic','yearly av. microbial respiration', &
                        '(kg-C/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'yrleach','yearly av. C leaching rate', &
                        '(kg-C/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'ayneetot','yearly TNEE','(kg-C/m**2/yr)',ndims, &
                        dimnames,istat)
            call inivar(idies,'ayanlit','yearly av. aboveground litter nitrogen', &
                        '(kg-N/m**2)',ndims,dimnames,istat)
            call inivar(idies,'aybnlit','yearly av. belowground litter nitrogen', &
                        '(kg-N/m**2)',ndims,dimnames,istat)
            call inivar(idies,'aynsoi','yearly av. soil nitrogen content', &
                        '(kg-N/m**2)',ndims,dimnames,istat)
            call inivar(idies,'ynleach','yearly av. N leaching rate', &
                        '(kg-N/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'ayprcp','yearly av. precipitation','(mm/yr)', &
                        ndims,dimnames,istat)
            call inivar(idies,'ayaet','yearly av. actual evapotranspiration', &
                        '(mm/yr)',ndims,dimnames,istat)
            call inivar(idies,'aytrans','yearly av. transpiration','(mm/yr)', &
                        ndims,dimnames,istat)
            call inivar(idies,'aysrunoff','yearly av. surface runoff', &
                        '(mm/yr)',ndims,dimnames,istat)
            call inivar(idies,'aydrainage','yearly av. drainage','(mm/yr)',ndims, &
                        dimnames,istat)
            call inivar(idies,'dwtot','yearly av. h20 change','(mm/yr)',ndims, &
                        dimnames,istat)
            call inivar(idies,'wtot','yearly av. h20 content','(mm/yr)',ndims, &
                        dimnames,istat)

            call inivar(idies,'aysolar','yearly av. incident solar radiation', &
                        '(w/m**2)',ndims,dimnames,istat)
            call inivar(idies,'ayreflect','yearly av. albedo','(fraction)',ndims, &
                        dimnames,istat)
            call inivar(idies,'ayirup','yearly av. upward ir radiation', &
                        '(w/m**2)',ndims,dimnames,istat)
            call inivar(idies,'ayirdown','yearly av. downward ir radiation', &
                        '(w/m**2)',ndims,dimnames,istat)
            call inivar(idies,'aysens','yearly av. sensible heat flux', &
                        '(w/m**2)',ndims,dimnames,istat)
            call inivar(idies,'aylatent','yearly av. latent heat flux', &
                        '(w/m**2)',ndims,dimnames,istat)
            call inivar(idies,'aywsoi','yearly av. 1m volumetric water content', &
                        '(fraction)',ndims,dimnames,istat)
            call inivar(idies,'aywisoi','yearly av. 1m volumetric ice content', &
                        '(fraction)',ndims,dimnames,istat)
            call inivar(idies,'aytsoi','yearly av. 1m soil temp.', &
                        '(fraction)',ndims,dimnames,istat)
            call inivar(idies,'ayvwc','yearly av. annual average 1m volumetric water content', &
                        '(fraction)',ndims,dimnames,istat)
            call inivar(idies,'ayawc','yearly av. annual average 1m plant available water content', &
                        '(fraction)',ndims,dimnames,istat)

            call inivar(idies,'firefac','yearly av. live root biomass', &
                        '(kg-C/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'ayrootbio','yearly av. live root biomass', &
                        '(kg-C/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'aynmintot','yearly av. total nitrogen mineralization', &
                        '(kg-N/m**2/yr)',ndims,dimnames,istat)
            call inivar(idies,'aycmic','yearly av. microbial biomass', &
                        '(kg-c/m**2)',ndims,dimnames,istat)
            call inivar(idies,'janxinwet','climatologic january precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'febxinwet','climatologic february precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'marxinwet','climatologic march precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'aprxinwet','climatologic april precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'mayxinwet','climatologic may precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'junxinwet','climatologic june precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'julxinwet','climatologic july precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'augxinwet','climatologic august precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'sepxinwet','climatologic september precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'octxinwet','climatologic october precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'novxinwet','climatologic november precipitation days', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'decxinwet','climatologic december precipitation days', &
                        '(days/month)',ndims,dimnames,istat)

! crops
         if(isimagro .gt. 0)then

            call inivar(idies,'gddpl15','accumulated growing degree days for 1 and half year for sugarcane', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'gddsgcp1','factor to get the hybrid planted for sugarcane plating, 1 and half year', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'gddsgcp2','factor to get the hybrid planted for sugarcane plating, 1 and half year', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'htmx1','maximum height attained by a crop during year', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'htmx2','maximum height attained by a crop during year', &
                        '(days/month)',ndims,dimnames,istat)

            call inivar(idies,'ik','counter for sugarcane crop one and half year GDD', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'cropy','sugarcane crop year since planting', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'cdays','days since planting in weather.f', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'ncyears','number of the ith crop year', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'ctot','total inorganic nitrogen in soil profile (kg n  m-2)', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'drntot','total inorganic nitrogen in soil profile (kg n  m-2)', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'cic3','intercellular co2 concentration - c3 crops (mol_co2/mol_air)', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'cic4','intercellular co2 concentration - c4 crops (mol_co2/mol_air)', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'csc3','leaf boundary layer co2 concentration - c3 crops (mol_co2/mol_air)', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'csc4','leaf boundary layer co2 concentration - c4 crops (mol_co2/mol_air)', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'ua',' wind speed (m s-1)', &
                        '(days/month)',ndims,dimnames,istat)

            call inivar(idies,'gddemerg','accumulated growing degree days at leaf emergence', &
                        '(days/month)',ndims,dimnames,istat)
            call inivar(idies,'a10tmin','10-day average minimum air temperature (K)', &
                        '(days/month)',ndims,dimnames,istat)

      endif


            call closefile(idies,istat)

         endif

         call openfile(idies,filen,istat)

      endif


!-----------------------------------------------------------------------
! writevar

! ===== xstore

      istart(3) = 1
      icount(3) = 3
      istart(4) = 1
      icount(4) = 1
      
      call writevar(filen,idies,'xstore',xstore,istart,icount,ftime,istat)

! ===== snow layers

      istart(3) = 1
      icount(3) = nsnolay
      istart(4) = 1
      icount(4) = mlpt
      istart(5) = 1
      icount(5) = 1
      
! temperature of snow layers
      call writevar(filen,idies,'tsno',tsno,istart,icount,ftime,istat)
! thickness of snow layers
      call writevar(filen,idies,'hsno',hsno,istart,icount,ftime,istat)
 

! ===== soil layers

      istart(3) = 1
      icount(3) = nsoilay
      istart(4) = 1
      icount(4) = mlpt
      istart(5) = 1
      icount(5) = 1

! temperature
      call writevar(filen,idies,'tsoi',tsoi,istart,icount,ftime,istat)
! ice content of soil
      call writevar(filen,idies,'wisoi',wisoi,istart,icount,ftime,istat)
! water content of soil
      call writevar(filen,idies,'wsoi',wsoi,istart,icount,ftime,istat)
!crops
     if(isimagro .gt. 0)then
      call writevar(filen,idies,'smsoil',smsoil,istart,icount,ftime,istat)
      call writevar(filen,idies,'smsoln',smsoln,istart,icount,ftime,istat)
     endif

! ===== Carbon pools and fluxes per pft

      istart(3) = 1
      icount(3) = npft
      istart(4) = 1
      icount(4) = mlpt
      istart(5) = 1
      icount(5) = 1
      
! carbon in leaf
      call writevar(filen,idies,'cbiol',cbiol,istart,icount,ftime,istat)
! carbon in wood
      call writevar(filen,idies,'cbiow',cbiow,istart,icount,ftime,istat)
! carbon in root
      call writevar(filen,idies,'cbior',cbior,istart,icount,ftime,istat)
!
      call writevar(filen,idies,'cbios',cbios,istart,icount,ftime,istat)
! current average npp per pft
      call writevar(filen,idies,'aynpp',aynpp,istart,icount,ftime,istat)
! current average gpp per pft
      call writevar(filen,idies,'aygpp',aygpp,istart,icount,ftime,istat)

!crops
   if(isimagro .gt. 0)then
     call writevar(filen,idies,'croplive',croplive,istart,icount,ftime,istat)
     call writevar(filen,idies,'gddplant',gddplant,istart,icount,ftime,istat)
     call writevar(filen,idies,'gddtsoi',gddtsoi,istart,icount,ftime,istat)
     call writevar(filen,idies,'idpp',idpp,istart,icount,ftime,istat)
     call writevar(filen,idies,'pstart',pstart,istart,icount,ftime,istat)
     call writevar(filen,idies,'cntops',cntops,istart,icount,ftime,istat)
     call writevar(filen,idies,'pdate',pdate,istart,icount,ftime,istat)
     call writevar(filen,idies,'hdate',hdate,istart,icount,ftime,istat)
     call writevar(filen,idies,'harvidx',harvidx,istart,icount,ftime,istat)
     call writevar(filen,idies,'croplaimx',croplaimx,istart,icount,ftime,istat)
     call writevar(filen,idies,'cropyld',cropyld,istart,icount,ftime,istat)
     call writevar(filen,idies,'dmyield',dmyield,istart,icount,ftime,istat)
     call writevar(filen,idies,'dmleaf',dmleaf,istart,icount,ftime,istat)
     call writevar(filen,idies,'dmstem',dmstem,istart,icount,ftime,istat)
     call writevar(filen,idies,'dmroot',dmroot,istart,icount,ftime,istat)
     call writevar(filen,idies,'dmresidue',dmresidue,istart,icount,ftime,istat)
     call writevar(filen,idies,'residuen',residuen,istart,icount,ftime,istat)
     call writevar(filen,idies,'nconcl',nconcl,istart,icount,ftime,istat)
     call writevar(filen,idies,'nconcs',nconcs,istart,icount,ftime,istat)
     call writevar(filen,idies,'nconcr',nconcr,istart,icount,ftime,istat)
     call writevar(filen,idies,'nconcg',nconcg,istart,icount,ftime,istat)
     call writevar(filen,idies,'cropn',cropn,istart,icount,ftime,istat)
     call writevar(filen,idies,'cropfixn',cropfixn,istart,icount,ftime,istat)
     call writevar(filen,idies,'cnrootvec',cnrootvec,istart,icount,ftime,istat)
     call writevar(filen,idies,'crmclim',crmclim,istart,icount,ftime,istat)
     call writevar(filen,idies,'crmact',crmact,istart,icount,ftime,istat)
     call writevar(filen,idies,'crmplant',crmplant,istart,icount,ftime,istat)
     call writevar(filen,idies,'aybprod',aybprod,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayabprod',ayabprod,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayrprod',ayrprod,istart,icount,ftime,istat)
     call writevar(filen,idies,'aylprod',aylprod,istart,icount,ftime,istat)
     call writevar(filen,idies,'aylprod',aylprod,istart,icount,ftime,istat)
     call writevar(filen,idies,'harvdate',harvdate,istart,icount,ftime,istat)
     call writevar(filen,idies,'cropplant',cropplant,istart,icount,ftime,istat)
     call writevar(filen,idies,'biomass',biomass,istart,icount,ftime,istat)
     call writevar(filen,idies,'hui',hui,istart,icount,ftime,istat)
     call writevar(filen,idies,'idop',idop,istart,icount,ftime,istat)
     call writevar(filen,idies,'fnleaf',fnleaf,istart,icount,ftime,istat)
     call writevar(filen,idies,'fnstem',fnstem,istart,icount,ftime,istat)
     call writevar(filen,idies,'stressn',stressn,istart,icount,ftime,istat)
     call writevar(filen,idies,'thrlai',thrlai,istart,icount,ftime,istat)
     call writevar(filen,idies,'grnfraccrop',grnfraccrop,istart,icount,ftime,istat)
     call writevar(filen,idies,'gddmaturity',gddmaturity,istart,icount,ftime,istat)
     call writevar(filen,idies,'adnpp',adnpp,istart,icount,ftime,istat)
     call writevar(filen,idies,'tnpp',tnpp,istart,icount,ftime,istat)
     call writevar(filen,idies,'plaimx',plaimx,istart,icount,ftime,istat)
     call writevar(filen,idies,'plai',plai,istart,icount,ftime,istat)
     call writevar(filen,idies,'totnuptake',totnuptake,istart,icount,ftime,istat)
     call writevar(filen,idies,'totnfix',totnfix,istart,icount,ftime,istat)
     call writevar(filen,idies,'cbiog',cbiog,istart,icount,ftime,istat)
     call writevar(filen,idies,'aroot',aroot,istart,icount,ftime,istat)
     call writevar(filen,idies,'aleaf',aleaf,istart,icount,ftime,istat)
     call writevar(filen,idies,'awood',awood,istart,icount,ftime,istat)
   endif

! ===== All other 2-d variables
      istart(3) = 1
      icount(3) = mlpt
      istart(4:) = 1
      icount(4:) = 1

! intrayear timestep index
      do n=lbeg,lend
         anytime(n) = nytime(1)
      enddo
      call writevar(filen,idies,'anytime',anytime,istart,icount,ftime,istat)
! fractional snow cover
! renamed fsnocov to fi which is the real var. name
      call writevar(filen,idies,'fi',fi,istart,icount,ftime,istat)
! sapwood fraction
      call writevar(filen,idies,'sapfrac',sapfrac,istart,icount,ftime,istat)
! leaf metabolic litter
      call writevar(filen,idies,'clitlm',clitlm,istart,icount,ftime,istat)
! leaf structural carbon litter
      call writevar(filen,idies,'clitls',clitls,istart,icount,ftime,istat)
! leaf lignin carbon litter
      call writevar(filen,idies,'clitll',clitll,istart,icount,ftime,istat)
! root metabolic litter
      call writevar(filen,idies,'clitrm',clitrm,istart,icount,ftime,istat)
! root structural litter
      call writevar(filen,idies,'clitrs',clitrs,istart,icount,ftime,istat)
! root lignin litter
      call writevar(filen,idies,'clitrl',clitrl,istart,icount,ftime,istat)
 ! woody metabolic litter
     call writevar(filen,idies,'clitwm',clitwm,istart,icount,ftime,istat)
! woody structural litter
     call writevar(filen,idies,'clitws',clitws,istart,icount,ftime,istat)
! woody lignin litter
     call writevar(filen,idies,'clitwl',clitwl,istart,icount,ftime,istat)
! annual leaf litterfall
     call writevar(filen,idies,'falll',falll,istart,icount,ftime,istat)
! annual fine root turnover 
     call writevar(filen,idies,'fallr',fallr,istart,icount,ftime,istat)
! annual wood turnover 
     call writevar(filen,idies,'fallw',fallw,istart,icount,ftime,istat)
! Vegetation type
     call writevar(filen,idies,'vegtype0',vegtype0,istart,icount,ftime,istat)
! Upper canopy fractional cover
     call writevar(filen,idies,'fu',fu,istart,icount,ftime,istat)
! Lower canopy fractional cover
     call writevar(filen,idies,'fl',fl,istart,icount,ftime,istat)
    if(isimagro .gt. 0)then
     call writevar(filen,idies,'gdd8',gdd8,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd10',gdd10,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd12',gdd12,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd0c',gdd0c,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd0clim',gdd0c,istart,icount,ftime,istat)
     call writevar(filen,idies,'ftot',ftot,istart,icount,ftime,istat)
   endif

!crops
   if(isimagro .gt. 0)then
     call writevar(filen,idies,'ctot',ctot,istart,icount,ftime,istat)
     call writevar(filen,idies,'drntot',drntot,istart,icount,ftime,istat)
     call writevar(filen,idies,'cic3',cic3,istart,icount,ftime,istat)
     call writevar(filen,idies,'cic4',cic4,istart,icount,ftime,istat)
     call writevar(filen,idies,'csc3',csc3,istart,icount,ftime,istat)
     call writevar(filen,idies,'csc4',csc4,istart,icount,ftime,istat)
     call writevar(filen,idies,'ua',ua,istart,icount,ftime,istat)
     call writevar(filen,idies,'dtnleach',dtnleach,istart,icount,ftime,istat)
   endif

! total microbial carbon
     call writevar(filen,idies,'totcmic',totcmic,istart,icount,ftime,istat)
! slow soil carbon, protected humus
     call writevar(filen,idies,'csoislop',csoislop,istart,icount,ftime,istat)
! slow soil carbon, nonprotected humus
     call writevar(filen,idies,'csoislon',csoislon,istart,icount,ftime,istat)
! passive soil carbon
     call writevar(filen,idies,'csoipas',csoipas,istart,icount,ftime,istat)
! average soil decomposition factor
     call writevar(filen,idies,'decomps',decomps,istart,icount,ftime,istat)
! average litter decomposition factor
     call writevar(filen,idies,'decompl',decompl,istart,icount,ftime,istat)
! growing degree days
! TODO: theoretically fixed vegetation runs does not require gdd0/5 gddthis0/5
!      AND pfts information to be written. Then we should disable writing this
!      on static vegetation runs!
     call writevar(filen,idies,'gdd0',gdd0,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd5',gdd5,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd0this',gdd0this,istart,icount,ftime,istat)
     call writevar(filen,idies,'gdd5this',gdd5this,istart,icount,ftime,istat)
! coldest monthly temperature
     call writevar(filen,idies,'tc',tc,istart,icount,ftime,istat)
     call writevar(filen,idies,'tcthis',tcthis,istart,icount,ftime,istat)
! warmest monthly temperature
     call writevar(filen,idies,'tw',tw,istart,icount,ftime,istat)
     call writevar(filen,idies,'twthis',twthis,istart,icount,ftime,istat)
! ice content of puddles
     call writevar(filen,idies,'wipud',wipud,istart,icount,ftime,istat)
! liquid content of puddles
     call writevar(filen,idies,'wpud',wpud,istart,icount,ftime,istat)
! intercepted water by lower canopy leaves and stems
     call writevar(filen,idies,'wliql',wliql,istart,icount,ftime,istat)
! intercepted water upper canopy stems 
     call writevar(filen,idies,'wliqs',wliqs,istart,icount,ftime,istat)
! intercepted water by upper canopy leaves
     call writevar(filen,idies,'wliqu',wliqu,istart,icount,ftime,istat)
! intercepted snow by lower canopy leaves and stems
     call writevar(filen,idies,'wsnol',wsnol,istart,icount,ftime,istat)
! intercepted snow by upper canopy stems
     call writevar(filen,idies,'wsnos',wsnos,istart,icount,ftime,istat)
! intercepted snow by upper canopy leaves
     call writevar(filen,idies,'wsnou',wsnou,istart,icount,ftime,istat)
! temperature of upper canopy leaves
     call writevar(filen,idies,'tu',tu,istart,icount,ftime,istat)
! temperature of upper canopy stems
     call writevar(filen,idies,'ts',ts,istart,icount,ftime,istat)
! temperature of lower canopy leaves
     call writevar(filen,idies,'tl',tl,istart,icount,ftime,istat)
! temperature at level z12
     call writevar(filen,idies,'t12',t12,istart,icount,ftime,istat)
! temperature at level z34
     call writevar(filen,idies,'t34',t34,istart,icount,ftime,istat)
! temperature of snow burried lower canopy
     call writevar(filen,idies,'tlsub',tlsub,istart,icount,ftime,istat)
! soil skin temperature
     call writevar(filen,idies,'tg',tg,istart,icount,ftime,istat)
! snow skin temperature
     call writevar(filen,idies,'ti',ti,istart,icount,ftime,istat)
! humidity at level z12
     call writevar(filen,idies,'q12',q12,istart,icount,ftime,istat)
! humidity at level z34
     call writevar(filen,idies,'q34',q34,istart,icount,ftime,istat)

! Physiology variables

! intercellular co2 concentration
     call writevar(filen,idies,'ciub',ciub,istart,icount,ftime,istat)
     call writevar(filen,idies,'ciuc',ciuc,istart,icount,ftime,istat)
     call writevar(filen,idies,'cils',cils,istart,icount,ftime,istat)
     call writevar(filen,idies,'cil3',cil3,istart,icount,ftime,istat)
     call writevar(filen,idies,'cil4',cil4,istart,icount,ftime,istat)
! leaf boundary layer co2 concentration
     call writevar(filen,idies,'csub',csub,istart,icount,ftime,istat)
     call writevar(filen,idies,'csuc',csuc,istart,icount,ftime,istat)
     call writevar(filen,idies,'csls',csls,istart,icount,ftime,istat)
     call writevar(filen,idies,'csl3',csl3,istart,icount,ftime,istat)
     call writevar(filen,idies,'csl4',csl4,istart,icount,ftime,istat)
! stomatal conductance
     call writevar(filen,idies,'gsub',gsub,istart,icount,ftime,istat)
     call writevar(filen,idies,'gsuc',gsuc,istart,icount,ftime,istat)
     call writevar(filen,idies,'gsls',gsls,istart,icount,ftime,istat)
     call writevar(filen,idies,'gsl3',gsl3,istart,icount,ftime,istat)
     call writevar(filen,idies,'gsl4',gsl4,istart,icount,ftime,istat)

! Phenology variables

! annual accumulated growing degree days for bud burst, upper canopy
     call writevar(filen,idies,'agddu',agddu,istart,icount,ftime,istat)
! annual accumulated growing degree days for bud burst, lower canopy
     call writevar(filen,idies,'agddl',agddl,istart,icount,ftime,istat)
! cold-phenology trigger for trees
     call writevar(filen,idies,'tempu',tempu,istart,icount,ftime,istat)
! cold-phenology trigger for grasses/shrubs
     call writevar(filen,idies,'templ',templ,istart,icount,ftime,istat)
! average daily air temperature
     call writevar(filen,idies,'td',td,istart,icount,ftime,istat)
! 10-day average daily air temperature
     call writevar(filen,idies,'a10td',a10td,istart,icount,ftime,istat)
! 10-day average canopy photosynthesis rate - broadleaf
     call writevar(filen,idies,'a10ancub',a10ancub,istart,icount,ftime,istat)
! 10-day average canopy photosynthesis rate - shrubs
     call writevar(filen,idies,'a10ancls',a10ancls,istart,icount,ftime,istat)
! 10-day average canopy photosynthesis rate - c4 grasses
     call writevar(filen,idies,'a10ancl4',a10ancl4,istart,icount,ftime,istat)
! 10-day average canopy photosynthesis rate - c3 grasses
     call writevar(filen,idies,'a10ancl3',a10ancl3,istart,icount,ftime,istat)
! 10-day average canopy scaling parameter - upper canopy
     call writevar(filen,idies,'a10scalparamu',a10scalparamu,istart,icount,ftime,istat)
! 10-day average canopy scaling parameter - lower canopy
     call writevar(filen,idies,'a10scalparaml',a10scalparaml,istart,icount,ftime,istat)
! 10-day average daylight - upper canopy
     call writevar(filen,idies,'a10daylightu',a10daylightu,istart,icount,ftime,istat)
! 10-day average daylight - lower canopy
     call writevar(filen,idies,'a10daylightl',a10daylightl,istart,icount,ftime,istat)
! drought-phenology trigger for trees
     call writevar(filen,idies,'dropu',dropu,istart,icount,ftime,istat)
! drought-phenology trigger for shrubs
     call writevar(filen,idies,'dropls',dropls,istart,icount,ftime,istat)
! drought-phenology trigger for c4 grasses
     call writevar(filen,idies,'dropl4',dropl4,istart,icount,ftime,istat)
! drought-phenology trigger for c3 grasses
     call writevar(filen,idies,'dropl3',dropl3,istart,icount,ftime,istat)
! LAI of upper canopy
     call writevar(filen,idies,'laiu',lai(:,2),istart,icount,ftime,istat)
! LAI of lower canopy
     call writevar(filen,idies,'lail',lai(:,1),istart,icount,ftime,istat)

! Summation variables

     call writevar(filen,idies,'ayalit',ayalit,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayblit',ayblit,istart,icount,ftime,istat)
     call writevar(filen,idies,'aycsoi',aycsoi,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayco2root',ayco2root,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayco2mic',ayco2mic,istart,icount,ftime,istat)
     call writevar(filen,idies,'yrleach',yrleach,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayneetot',ayneetot,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayanlit',ayanlit,istart,icount,ftime,istat)
     call writevar(filen,idies,'aybnlit',aybnlit,istart,icount,ftime,istat)
     call writevar(filen,idies,'aynsoi',aynsoi,istart,icount,ftime,istat)
     call writevar(filen,idies,'ynleach',ynleach,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayprcp',ayprcp,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayaet',ayaet,istart,icount,ftime,istat)
     call writevar(filen,idies,'aytrans',aytrans,istart,icount,ftime,istat)
     call writevar(filen,idies,'aysrunoff',aysrunoff,istart,icount,ftime,istat)
     call writevar(filen,idies,'aydrainage',aydrainage,istart,icount,ftime,istat)
     call writevar(filen,idies,'dwtot',dwtot,istart,icount,ftime,istat)
     call writevar(filen,idies,'wtot',wtot,istart,icount,ftime,istat)
     call writevar(filen,idies,'aysolar',aysolar,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayreflect',ayreflect,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayirup',ayirup,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayirdown',ayirdown,istart,icount,ftime,istat)
     call writevar(filen,idies,'aysens',aysens,istart,icount,ftime,istat)
     call writevar(filen,idies,'aylatent',aylatent,istart,icount,ftime,istat)
     call writevar(filen,idies,'aywsoi',aywsoi,istart,icount,ftime,istat)
     call writevar(filen,idies,'aywisoi',aywisoi,istart,icount,ftime,istat)
     call writevar(filen,idies,'aytsoi',aytsoi,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayvwc',ayvwc,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayawc',ayawc,istart,icount,ftime,istat)
     call writevar(filen,idies,'firefac',firefac,istart,icount,ftime,istat)
     call writevar(filen,idies,'ayrootbio',ayrootbio,istart,icount,ftime,istat)
     call writevar(filen,idies,'aynmintot',aynmintot,istart,icount,ftime,istat)
     call writevar(filen,idies,'aycmic',aycmic,istart,icount,ftime,istat)
!xinwet
     call writevar(filen,idies,'janxinwet',xinwet(:,1),istart,icount,ftime,istat)
     call writevar(filen,idies,'febxinwet',xinwet(:,2),istart,icount,ftime,istat)
     call writevar(filen,idies,'marxinwet',xinwet(:,3),istart,icount,ftime,istat)
     call writevar(filen,idies,'aprxinwet',xinwet(:,4),istart,icount,ftime,istat)
     call writevar(filen,idies,'mayxinwet',xinwet(:,5),istart,icount,ftime,istat)
     call writevar(filen,idies,'junxinwet',xinwet(:,6),istart,icount,ftime,istat)
     call writevar(filen,idies,'julxinwet',xinwet(:,7),istart,icount,ftime,istat)
     call writevar(filen,idies,'augxinwet',xinwet(:,8),istart,icount,ftime,istat)
     call writevar(filen,idies,'sepxinwet',xinwet(:,9),istart,icount,ftime,istat)
     call writevar(filen,idies,'octxinwet',xinwet(:,10),istart,icount,ftime,istat)
     call writevar(filen,idies,'novxinwet',xinwet(:,11),istart,icount,ftime,istat)
     call writevar(filen,idies,'decxinwet',xinwet(:,12),istart,icount,ftime,istat)

!crops
   if(isimagro .gt. 0)then
     call writevar(filen,idies,'gddpl15',gddpl15,istart,icount,ftime,istat)
     call writevar(filen,idies,'gddsgcp1',gddsgcp(:,1),istart,icount,ftime,istat)
     call writevar(filen,idies,'gddsgcp2',gddsgcp(:,2),istart,icount,ftime,istat)
     call writevar(filen,idies,'htmx1',htmx(:,1),istart,icount,ftime,istat)
     call writevar(filen,idies,'htmx2',htmx(:,2),istart,icount,ftime,istat)
     call writevar(filen,idies,'ik',ik,istart,icount,ftime,istat)
     call writevar(filen,idies,'cropy',cropy,istart,icount,ftime,istat)
     call writevar(filen,idies,'cdays',cdays,istart,icount,ftime,istat)
     call writevar(filen,idies,'ncyears',ncyears,istart,icount,ftime,istat)
     call writevar(filen,idies,'aplantn',aplantn,istart,icount,ftime,istat)
     call writevar(filen,idies,'gddemerg',gddemerg,istart,icount,ftime,istat)
     call writevar(filen,idies,'a10tmin',a10tmin,istart,icount,ftime,istat)
   endif

      istart(3) = 1
      icount(3) = (ecpft-scpft+1)
      istart(ndims+1) = mstep
      icount(ndims+1) = 1
   if(isimagro .gt. 0)then
      call writevar(filen,idies,'cropout(:,4,35)',cropout(:,4:,35),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,5,35)',cropout(:,5:,35),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,6,35)',cropout(:,6:,35),istart,icount,ftime,istat)

      call writevar(filen,idies,'cropout(:,13:,50)',cropout(:,13:,50),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,14:,50)',cropout(:,14:,50),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,15:,50)',cropout(:,15:,50),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,16:,50)',cropout(:,16:,50),istart,icount,ftime,istat)


      call writevar(filen,idies,'cropout(:,13:,52)',cropout(:,13:,52),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,14:,52)',cropout(:,14:,52),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,15:,52)',cropout(:,15:,52),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,16:,52)',cropout(:,16:,52),istart,icount,ftime,istat)

      call writevar(filen,idies,'cropout(:,13:,51)',cropout(:,13:,51),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,14:,51)',cropout(:,14:,51),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,15:,51)',cropout(:,15:,51),istart,icount,ftime,istat)
      call writevar(filen,idies,'cropout(:,16:,51)',cropout(:,16:,51),istart,icount,ftime,istat)
  endif
     ! close file
     call closefile(idies,istat)

     ! write control file
     if (myid .eq. 0) then
        open(unit=8, file=trim(outdir)//'/inland-restart_control.dat')
        write(8,*) iyear, imonthp
        close(8)
     endif
     
     ! restore env_floatout
     env_floatout = tmp_floatout
     
     return
end subroutine wrestart
