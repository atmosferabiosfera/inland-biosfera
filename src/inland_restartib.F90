#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! ---------------------------------------------------------------------
subroutine restartib (iyearp, imonthp)
! ---------------------------------------------------------------------
! reads in restart files, initializes some variables
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_control
      use inland_comveg
      use inland_comsoi
      use inland_comsno
      use inland_combcs
      use inland_comwork
      use inland_comsum
      use inland_comatm, only: xstore
      use inland_comnitr
      use inland_comcrop
      use inland_com1d, only:straml, strahl, richl, richu, stramu, strahu, &
                        raing, rainl, rainu, pfluxu, pfluxs, pfluxl, za

      implicit none
! ------------------------------Arguments--------------------------------
! Input arguments
      integer :: iyearp, imonthp

! this subroutine reads the restart values of:
!
! fi      = fractional snow cover
! tsno    = temperature of snow
! hsno    = snow depth
! tsoi    = soil temperature
! wisoi   = soil ice content
! wsoi    = soil moisture content
! cbiol   = carbon in leaf biomass pool
! cbiow   = carbon in woody biomass pool
! cbior   = carbon in fine root biomass pool
! sapfrac = sapwood fraction
! clitlm  = leaf metabolic litter
! clitls  = leaf structural litter
! clitll  = leaf lignin litter
! clitrm  = root metabolic litter
! clitrs  = root structural litter
! clitrl  = root lignin litter
! clitwm  = woody metabolic litter
! clitws  = woody structural litter
! clitwl  = woody lignin litter
! falll   = annual leaf litterfall
! fallr   = annual fine root turnover
! fallw   = annual woody turnover
! vegtype0= vegetation type
! fu      = fraction of overall area covered by upper canopy
! fl      = fraction of snow-free area covered by lower  canopy
! totcmic = total microbial carbon
! csoislop= slow soil carbon, protected humus
! csoislon= slow soil carbon, nonprotected humus
! csoipas = passive soil carbon
! decomps = soil decomposition factor
! decompl = litter decomposition factor
! gdd0    = growing degree days 0
! gdd5    = growing degree days 5
! tc      = coldest monthly temperature
! tw      = warmest monthly temperature
! gdd0this= growing degree days 0 of the current year
! gdd5this= growing degree days 5 of the current year
! tcthis  = coldest monthly temperature of the current year
! twthis  = warmest monthly temperature of the current year
! wipud   = ice content of puddles per soil area
! wpud    = liquid content of puddles per soil area
! wliql   = intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
! wliqs   = intercepted liquid h2o on upper canopy stem area (kg m-2)
! wliqu   = intercepted liquid h2o on upper canopy leaf area (kg m-2)
! wsnol   = intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
! wsnos   = intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
! wsnou   = intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
! tu      = temperature of upper canopy leaves (K)
! ts      = temperature of upper canopy stems (K)
! tl      = temperature of lower canopy leaves & stems(K)
! t12     = air temperature at z12 (K)
! t34     = air temperature at z34 (K)
! tlsub   = temperature of lower canopy vegetation buried by snow (K)
! tg      = soil skin temperature (K)
! ti      = snow skin temperature (K)
! q12     = specific humidity of air at z12
! q34     = specific humidity of air at z34
!
! Physiology
! ciub    = intercellular co2 concentration - broadleaf (mol_co2/mol_air)
! ciuc    = intercellular co2 concentration - conifer   (mol_co2/mol_air)
! cils    = intercellular co2 concentration - shrubs    (mol_co2/mol_air)
! cil3    = intercellular co2 concentration - c3 plants (mol_co2/mol_air)
! cil4    = intercellular co2 concentration - c4 plants (mol_co2/mol_air)
! csub    = leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
! csuc    = leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
! csls    = leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
! csl3    = leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
! csl4    = leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
! gsub    = upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
! gsuc    = upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
! gsls    = lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
! gsl3    = lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
! gsl4    = lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
!
! Phenology :
! agddu   = annual accumulated growing degree days for bud burst, upper canopy
! agddl   = annual accumulated growing degree days for bud burst, lower canopy
! tempu   = cold-phenology trigger for trees
! templ   = cold-phenology trigger for grasses/shrubs
! td      = avg daily temp
! a10td    = 10-day avg daily temp
! a10ancub = 10-day average canopy photosynthesis rate - broadleaf
! a10ancuc = 10-day average canopy photosynthesis rate - conifer
! a10ancls = 10-day average canopy photosynthesis rate - shrubs
! a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
! a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
! a10scalparamu = 10-day average canopy scaling parameter - upper canopy
! a10scalparaml = 10-day average canopy scaling parameter - lower canopy
! a10daylightu = 10-day average daylight - upper canopy
! a10daylightl = 10-day average daylight - lower canopy
! dropu   = drought-phenology trigger for trees
! dropls  = drought-phenology trigger for shrubs
! dropl4  = drought-phenology trigger for c4 grasses
! dropl3  = drought-phenology trigger for c3 grasses
! lail = lai(i,1) : lower canopy lai
! laiu = lai(i,2) : upper canopy lai
! (NOTE: a10ancuc is not used at this point, so its restart entry
! is commented out)
!
! Summation variables used in dynaveg, gdiag and vdiag
!
! aynpp = annual total npp per pft (kg-c/m**2/yr)
! aygpp = annual total gpp per pft (kg-c/m**2/yr)
! ayalit = yearly av. aboveground litter (kg-c/m**2)
! ayblit = yearly av. belowground litter (kg-c/m**2)
! aycsoi = yearly av. soil carbon content (kg-c/m**2)
! ayco2root = yearly av. root respiration (kg-C/m**2/yr)
! ayco2mic = yearly av. microbial respiration (kg-C/m**2/yr)
! yrleach = yearly av. C leaching rate (kg_C m-2/yr)
! ayanlit = yearly av. aboveground litter nitrogen (kg-N/m**2)
! aybnlit = yearly av. belowground litter nitrogen (kg-N/m**2)
! aynsoi = yearly av. total soil nitrogen (kg-N/m**2)
! ynleach = yearly av. N leaching rate (kg_N m-2/yr)
! ayprcp = yearly av. precipitation (mm/year)
! ayaet  = yearly av. aet (mm/yr)
! aytrans = yearly av. transpiration (mm/yr)
! aysrunoff = yearly av. surface runoff (mm/yr)
! aydrainage = yearly av. drainage (mm/yr)
! dwtot = yearly av. change in total water content (mm/yr)
! wtot = yearly h20 content (mm/yr)
! aysolar = annual average incident solar radiation (w/m**2)
! ayreflect = annual average reflected radiation  (W/m**2)
! ayirup = annual average upward ir radiation (w/m**2)
! ayirdown = annual average downward ir radiation (w/m**2)
! aysens = annual average sensible heat flux (w/m**2)
! aylatent = annual average latent heat flux (w/m**2)
!
! aywsoi = yearly av. soil moisture content
! aywisoi = yearly av. soil ice content
! aytsoi = yearly av.  temperature
! ayvwc = annual average 1m volumetric water content (fraction)
! ayawc = annual average 1m plant-available water content (fraction)
! firefac = annual average fire factor
! ayrootbio = annual average live root biomass (kg-C/m**2/yr)
! aynmintot = annual total nitrogen mineralization (kg-N/m**2/yr)
! aycmic = annual av. microbial biomass (kg-c/m**2)
!
! local variables
      integer :: istart(5), icount(5), istat      ! for reading restart vars
      integer :: i, nlevel, ilpt
      integer :: lenchrib, lf, nday

      character*20 fdir  ! used to construct odd/even file names
      character*1024 filen ! file name
      character*4  chyear
      character*2  chmonth
      real*8 ftime         ! real form of nday
      integer ndims        ! number of dims (including tile, not including pft)

      character*80 varname

      data istart / 1,1,1,1,1 /, icount / nlon,nlat,1,1,1 /
! ---------------------------------------------------------------------
      icount(1) = nlonsub
      icount(2) = nlatsub

      chyear = '0000'
      if (iyearp .gt. 1000) then
         write(chyear(1:4),'(i4)') iyearp
      else if (iyearp .lt. 10) then
         write(chyear(4:4),'(i1)') iyearp
      else if (iyearp .lt. 100) then
         write(chyear(3:4),'(i2)') iyearp
      else
         write(chyear(2:4),'(i3)') iyearp
      endif

      chmonth = '00'
      if (imonthp .gt. 9) then
         write(chmonth(1:2),'(i2)') imonthp
      else
         write(chmonth(2:2),'(i1)') imonthp
      endif

      fdir = 'restart/'
      lf = lenchrib(fdir)

      ftime = nday - ndaypy + 1

! how many dims?
      if ( mlpt .gt. 1 ) then
         ndims = 4
      else
         ndims = 3
      end if

!-----------------------------------------------------------------------
! readvar

! by Etienne on 06/01/2013
! one restart file with all vars and dims!
! TODO what about if (myid .eq. 0) ???

      filen = trim(outdir)//'/inland-restart-'//chyear//'_'//chmonth//'.nc'
      
! ===== xstore

      ! loop for all xstorelevel - write each var, one layer at a time
      ! TODO investigate if this causes performance problems
       do nlevel = 1,3

          istart(3) = nlevel
          icount(3) = 1
          call readvar(filen,'xstore',xstore(:,nlevel),istart,icount,1,istat)

       end do ! xstore


! ===== snow layers

      ! loop for all nsnolay - write each var, one layer at a time
      ! TODO investigate if this causes performance problems

       do nlevel = 1,nsnolay

          istart(3) = nlevel
          icount(3) = 1
          istart(4) = 1
          icount(4) = mlpt

         call readvar(filen,'tsno',tsno(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'hsno',hsno(:,nlevel),istart,icount,mlpt,istat)
      end do


! ===== soil layers

      ! loop for all nsoilay - write each var, one layer at a time
      ! TODO investigate if this causes performance problems
      do nlevel = 1,nsoilay

         istart(3) = nlevel
         icount(3) = 1
         istart(4) = 1
         icount(4) = mlpt

         call readvar(filen,'tsoi',tsoi(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'wisoi',wisoi(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'wsoi',wsoi(:,nlevel),istart,icount,mlpt,istat)
!crops
       if(isimagro .gt. 0) then
         call readvar(filen,'smsoil',smsoil(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'smsoln',smsoln(:,nlevel),istart,icount,mlpt,istat)
       endif
      end do ! nsoilay
 


! ===== Carbon pools and fluxes per pft

      ! loop for all pfts - read each var, one pft at a time
      ! TODO investigate if this causes performance problems
     do nlevel = 1,npft

         istart(3) = nlevel
         icount(3) = 1
         istart(4) = 1
         icount(4) = mlpt

         call readvar(filen,'cbiol',cbiol(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cbiow',cbiow(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cbior',cbior(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cbios',cbios(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'aynpp',aynpp(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'aygpp',aygpp(:,nlevel),istart,icount,mlpt,istat)

! crops
       if(isimagro .gt. 0) then
         call readvar(filen,'croplive',croplive(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'gddplant',gddplant(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'gddtsoi',gddtsoi(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'idpp',idpp(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'pstart',pstart(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cntops',cntops(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'pdate',pdate(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'hdate',hdate(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'harvidx',harvidx(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'croplaimx',croplaimx(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cropyld',cropyld(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'dmyield',dmyield(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'dmleaf',dmleaf(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'dmstem',dmstem(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'dmroot',dmroot(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'dmresidue',dmresidue(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'residuen',residuen(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'nconcl',nconcl(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'nconcs',nconcs(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'nconcr',nconcr(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'nconcg',nconcg(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cropn',cropn(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cnrootvec',cnrootvec(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'crmclim',crmclim(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'crmact',crmact(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'crmplant',crmplant(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'aybprod',aybprod(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'ayabprod',ayabprod(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'ayrprod',ayrprod(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'aylprod',aylprod(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'harvdate',harvdate(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cropplant',cropplant(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'biomass',biomass(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'hui',hui(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'idop',idop(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'fnleaf',fnleaf(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'fnstem',fnstem(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'stressn',stressn(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'thrlai',thrlai(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'grnfraccrop',grnfraccrop(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'gddmaturity',gddmaturity(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'adnpp',adnpp(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'tnpp',tnpp(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'plaimx',plaimx(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'plai',plai(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'totnuptake',totnuptake(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'totnfix',totnfix(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'cbiog',cbiog(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'aroot',aroot(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'aleaf',aleaf(:,nlevel),istart,icount,mlpt,istat)
         call readvar(filen,'awood',awood(:,nlevel),istart,icount,mlpt,istat)
       endif
      end do ! npft

! ===== All other 2-d variables
      istart(3) = 1
      icount(3) = mlpt
      istart(4:) = 1
      icount(4:) = 1

! intrayear timestep index
! TODO test this
      varname = 'anytime'
      call readvar(filen,'anytime',anytime,istart,icount,0,istat)

      do i=1,numlv
         nytime(i) = int( anytime(lbeg) + 0.1 )
      end do

      ! renamed fsnocov to fi which is the real var. name
      call readvar(filen,'fi',fi,istart,icount,mlpt,istat)
      call readvar(filen,'sapfrac',sapfrac,istart,icount,mlpt,istat)
      call readvar(filen,'clitlm',clitlm,istart,icount,mlpt,istat)
      call readvar(filen,'clitls',clitls,istart,icount,mlpt,istat)
      call readvar(filen,'clitll',clitll,istart,icount,mlpt,istat)
      call readvar(filen,'clitrm',clitrm,istart,icount,mlpt,istat)
      call readvar(filen,'clitrs',clitrs,istart,icount,mlpt,istat)
      call readvar(filen,'clitrl',clitrl,istart,icount,mlpt,istat)
      call readvar(filen,'clitwm',clitwm,istart,icount,mlpt,istat)
      call readvar(filen,'clitws',clitws,istart,icount,mlpt,istat)
      call readvar(filen,'clitwl',clitwl,istart,icount,mlpt,istat)
      call readvar(filen,'falll',falll,istart,icount,mlpt,istat)
      call readvar(filen,'fallr',fallr,istart,icount,mlpt,istat)
      call readvar(filen,'fallw',fallw,istart,icount,mlpt,istat)
      call readvar(filen,'vegtype0',vegtype0,istart,icount,mlpt,istat)
      call readvar(filen,'fu',fu,istart,icount,mlpt,istat)
      call readvar(filen,'fl',fl,istart,icount,mlpt,istat)
      call readvar(filen,'totcmic',totcmic,istart,icount,mlpt,istat)
      call readvar(filen,'csoislop',csoislop,istart,icount,mlpt,istat)
      call readvar(filen,'csoislon',csoislon,istart,icount,mlpt,istat)
      call readvar(filen,'csoipas',csoipas,istart,icount,mlpt,istat)
      call readvar(filen,'decomps',decomps,istart,icount,mlpt,istat)
      call readvar(filen,'decompl',decompl,istart,icount,mlpt,istat)
!crops
    if(isimagro .gt. 0) then
      call readvar(filen,'gdd8',gdd8,istart,icount,ftime,istat)
      call readvar(filen,'gdd10',gdd10,istart,icount,ftime,istat)
      call readvar(filen,'gdd12',gdd12,istart,icount,ftime,istat)
      call readvar(filen,'gdd0c',gdd0c,istart,icount,ftime,istat)
      call readvar(filen,'gdd0clim',gdd0,istart,icount,ftime,istat)
      call readvar(filen,'ftot',ftot,istart,icount,ftime,istat)
      call readvar(filen,'ctot',ctot,istart,icount,mlpt,istat)
      call readvar(filen,'drntot',drntot,istart,icount,mlpt,istat)
      call readvar(filen,'cic3',cic3,istart,icount,mlpt,istat)
      call readvar(filen,'cic4',cic4,istart,icount,mlpt,istat)
      call readvar(filen,'csc3',csc3,istart,icount,mlpt,istat)
      call readvar(filen,'csc4',csc4,istart,icount,mlpt,istat)
      call readvar(filen,'ua',ua,istart,icount,mlpt,istat)
      call readvar(filen,'dtnleach',dtnleach,istart,icount,mlpt,istat)
    endif

!   read only in the case of dynamic vegetation restart or fixed vegetation 
! using the ccm3 climate
!   In the other cases,  climanl in inland_main_offline calculates gdd0, gdd5, tc, tw, 
! tcmin from climatology
! TODO what about tcmin???

!   The offline model always calls climanl once regardless of dynaveg or restart
! - fzm
!#ifdef COUPLED
!      if ((isimveg .gt. 0) .or. (ccmexist .gt. 0)) then
!   Modified because of read in restartib subroutine 
      if (isimveg .ge. 1) then

         call readvar(filen,'gdd0',gdd0,istart,icount,mlpt,istat)
         call readvar(filen,'gdd5',gdd5,istart,icount,mlpt,istat)
         call readvar(filen,'tc',tc,istart,icount,mlpt,istat)
         call readvar(filen,'tw',tw,istart,icount,mlpt,istat)

      end if ! isimveg
!#endif /* COUPLED */

      call readvar(filen,'gdd0this',gdd0this,istart,icount,mlpt,istat)
      call readvar(filen,'gdd5this',gdd5this,istart,icount,mlpt,istat)
      call readvar(filen,'tcthis',tcthis,istart,icount,mlpt,istat) 
      call readvar(filen,'twthis',twthis,istart,icount,mlpt,istat)
      call readvar(filen,'wipud',wipud,istart,icount,mlpt,istat)
      call readvar(filen,'wpud',wpud,istart,icount,mlpt,istat)
      call readvar(filen,'wliql',wliql,istart,icount,mlpt,istat)
      call readvar(filen,'wliqs',wliqs,istart,icount,mlpt,istat)
      call readvar(filen,'wliqu',wliqu,istart,icount,mlpt,istat)
      call readvar(filen,'wsnol',wsnol,istart,icount,mlpt,istat)
      call readvar(filen,'wsnos',wsnos,istart,icount,mlpt,istat)
      call readvar(filen,'wsnou',wsnou,istart,icount,mlpt,istat)
      call readvar(filen,'tu',tu,istart,icount,mlpt,istat)
      call readvar(filen,'ts',ts,istart,icount,mlpt,istat)
      call readvar(filen,'tl',tl,istart,icount,mlpt,istat)
      call readvar(filen,'t12',t12,istart,icount,mlpt,istat)
      call readvar(filen,'t34',t34,istart,icount,mlpt,istat)
      call readvar(filen,'tg',tg,istart,icount,mlpt,istat)
      call readvar(filen,'ti',ti,istart,icount,mlpt,istat)
      call readvar(filen,'tlsub',tlsub,istart,icount,mlpt,istat)
      call readvar(filen,'q12',q12,istart,icount,mlpt,istat)
      call readvar(filen,'q34',q34,istart,icount,mlpt,istat)

! Physiology variables

      call readvar(filen,'ciub',ciub,istart,icount,mlpt,istat)
      call readvar(filen,'ciuc',ciuc,istart,icount,mlpt,istat)
      call readvar(filen,'cils',cils,istart,icount,mlpt,istat)
      call readvar(filen,'cil3',cil3,istart,icount,mlpt,istat)
      call readvar(filen,'cil4',cil4,istart,icount,mlpt,istat)
      call readvar(filen,'csub',csub,istart,icount,mlpt,istat)
      call readvar(filen,'csuc',csuc,istart,icount,mlpt,istat)
      call readvar(filen,'csls',csls,istart,icount,mlpt,istat)
      call readvar(filen,'csl3',csl3,istart,icount,mlpt,istat)
      call readvar(filen,'csl4',csl4,istart,icount,mlpt,istat)
      call readvar(filen,'gsub',gsub,istart,icount,mlpt,istat)
      call readvar(filen,'gsuc',gsuc,istart,icount,mlpt,istat)
      call readvar(filen,'gsls',gsls,istart,icount,mlpt,istat)
      call readvar(filen,'gsl3',gsl3,istart,icount,mlpt,istat)
      call readvar(filen,'gsl4',gsl4,istart,icount,mlpt,istat)

! Phenology variables

      call readvar(filen,'agddu',agddu,istart,icount,mlpt,istat)
      call readvar(filen,'agddl',agddl,istart,icount,mlpt,istat)
      call readvar(filen,'tempu',tempu,istart,icount,mlpt,istat)
      call readvar(filen,'templ',templ,istart,icount,mlpt,istat)
      call readvar(filen,'td',td,istart,icount,mlpt,istat)
      call readvar(filen,'a10td',a10td,istart,icount,mlpt,istat)
      call readvar(filen,'a10ancub',a10ancub,istart,icount,mlpt,istat)
      call readvar(filen,'a10ancls',a10ancls,istart,icount,mlpt,istat)
      call readvar(filen,'a10ancl4',a10ancl4,istart,icount,mlpt,istat)
      call readvar(filen,'a10ancl3',a10ancl3,istart,icount,mlpt,istat)
      call readvar(filen,'a10scalparamu',a10scalparamu,istart,icount,mlpt,istat)
      call readvar(filen,'a10scalparaml',a10scalparaml,istart,icount,mlpt,istat)
      call readvar(filen,'a10daylightu',a10daylightu,istart,icount,mlpt,istat)
      call readvar(filen,'a10daylightl',a10daylightl,istart,icount,mlpt,istat)
      call readvar(filen,'dropu',dropu,istart,icount,mlpt,istat)
      call readvar(filen,'dropls',dropls,istart,icount,mlpt,istat)
      call readvar(filen,'dropl4',dropl4,istart,icount,mlpt,istat)
      call readvar(filen,'dropl3',dropl3,istart,icount,mlpt,istat)
      call readvar(filen,'laiu',lai(:,2),istart,icount,mlpt,istat)
      call readvar(filen,'lail',lai(:,1),istart,icount,mlpt,istat)

 ! Summation variables

      call readvar(filen,'ayalit',ayalit,istart,icount,mlpt,istat)
      call readvar(filen,'ayblit',ayblit,istart,icount,mlpt,istat)
      call readvar(filen,'aycsoi',aycsoi,istart,icount,mlpt,istat)
      call readvar(filen,'ayco2root',ayco2root,istart,icount,mlpt,istat)
      call readvar(filen,'ayco2mic',ayco2mic,istart,icount,mlpt,istat)
      call readvar(filen,'yrleach',yrleach,istart,icount,mlpt,istat)
      call readvar(filen,'ayneetot',ayneetot,istart,icount,mlpt,istat)
      call readvar(filen,'ayanlit',ayanlit,istart,icount,mlpt,istat)
      call readvar(filen,'aybnlit',aybnlit,istart,icount,mlpt,istat)
      call readvar(filen,'aynsoi',aynsoi,istart,icount,mlpt,istat)
      call readvar(filen,'ynleach',ynleach,istart,icount,mlpt,istat)
      call readvar(filen,'ayprcp',ayprcp,istart,icount,mlpt,istat)
      call readvar(filen,'ayaet',ayaet,istart,icount,mlpt,istat)
      call readvar(filen,'aytrans',aytrans,istart,icount,mlpt,istat)
      call readvar(filen,'aysrunoff',aysrunoff,istart,icount,mlpt,istat)
      call readvar(filen,'aydrainage',aydrainage,istart,icount,mlpt,istat)
      call readvar(filen,'dwtot',dwtot,istart,icount,mlpt,istat)
      call readvar(filen,'wtot',wtot,istart,icount,mlpt,istat)
      call readvar(filen,'aysolar',aysolar,istart,icount,mlpt,istat)
      call readvar(filen,'ayreflect',ayreflect,istart,icount,mlpt,istat)
      call readvar(filen,'ayirup',ayirup,istart,icount,mlpt,istat)
      call readvar(filen,'ayirdown',ayirdown,istart,icount,mlpt,istat)
      call readvar(filen,'aysens',aysens,istart,icount,mlpt,istat)
      call readvar(filen,'aylatent',aylatent,istart,icount,mlpt,istat)
      call readvar(filen,'aywsoi',aywsoi,istart,icount,mlpt,istat)
      call readvar(filen,'aywisoi',aywisoi,istart,icount,mlpt,istat)
      call readvar(filen,'aytsoi',aytsoi,istart,icount,mlpt,istat)
      call readvar(filen,'ayvwc',ayvwc,istart,icount,mlpt,istat)
      call readvar(filen,'ayawc',ayawc,istart,icount,mlpt,istat)
      call readvar(filen,'firefac',firefac,istart,icount,mlpt,istat)
      call readvar(filen,'ayrootbio',ayrootbio,istart,icount,mlpt,istat)
      call readvar(filen,'aynmintot',aynmintot,istart,icount,mlpt,istat)
      call readvar(filen,'aycmic',aycmic,istart,icount,mlpt,istat)
!xinwet
      call readvar(filen,'janxinwet',xinwet(:,1),istart,icount,mlpt,istat)
      call readvar(filen,'febxinwet',xinwet(:,2),istart,icount,mlpt,istat)
      call readvar(filen,'marxinwet',xinwet(:,3),istart,icount,mlpt,istat)
      call readvar(filen,'aprxinwet',xinwet(:,4),istart,icount,mlpt,istat)
      call readvar(filen,'mayxinwet',xinwet(:,5),istart,icount,mlpt,istat)
      call readvar(filen,'junxinwet',xinwet(:,6),istart,icount,mlpt,istat)
      call readvar(filen,'julxinwet',xinwet(:,7),istart,icount,mlpt,istat)
      call readvar(filen,'augxinwet',xinwet(:,8),istart,icount,mlpt,istat)
      call readvar(filen,'sepxinwet',xinwet(:,9),istart,icount,mlpt,istat)
      call readvar(filen,'octxinwet',xinwet(:,10),istart,icount,mlpt,istat)
      call readvar(filen,'novxinwet',xinwet(:,11),istart,icount,mlpt,istat)
      call readvar(filen,'decxinwet',xinwet(:,12),istart,icount,mlpt,istat)

! crops
    if(isimagro .gt. 0) then
      call readvar(filen,'ud',ud,istart,icount,mlpt,istat)
      call readvar(filen,'za',za,istart,icount,mlpt,istat)
      call readvar(filen,'consdays',consdays,istart,icount,mlpt,istat)
      call readvar(filen,'gddpl15',gddpl15,istart,icount,mlpt,istat)
      call readvar(filen,'gddsgcp1',gddsgcp(:,1),istart,icount,mlpt,istat)
      call readvar(filen,'gddsgcp2',gddsgcp(:,2),istart,icount,mlpt,istat)
      call readvar(filen,'htmx1',htmx(:,1),istart,icount,mlpt,istat)
      call readvar(filen,'htmx2',htmx(:,2),istart,icount,mlpt,istat)
      call readvar(filen,'ik',ik,istart,icount,mlpt,istat)
      call readvar(filen,'cropy',cropy,istart,icount,mlpt,istat)
      call readvar(filen,'cdays',cdays,istart,icount,mlpt,istat)
      call readvar(filen,'ncyears',ncyears,istart,icount,mlpt,istat)
      call readvar(filen,'aplantn',aplantn,istart,icount,mlpt,istat)
      call readvar(filen,'gddemerg',gddemerg,istart,icount,mlpt,istat)
      call readvar(filen,'a10tmin',a10tmin,istart,icount,mlpt,istat)
    endif

      istart(3) = 1
      icount(3) = (ecpft-scpft+1)
      istart(ndims+1) = 1
      icount(ndims+1) = 1

       if(isimagro .gt. 0) then
         call readvar(filen,'cropout(:,4,35)',cropout(:,4,35),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,5,35)',cropout(:,5,35),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,6,35)',cropout(:,6,35),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,13:,50)',cropout(:,13:,50),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,14:,50)',cropout(:,14:,50),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,15:,50)',cropout(:,15:,50),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,16:,50)',cropout(:,16:,50),istart,icount,mlpt,istat)

         call readvar(filen,'cropout(:,13:,52)',cropout(:,13:,52),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,14:,52)',cropout(:,14:,52),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,15:,52)',cropout(:,15:,52),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,16:,52)',cropout(:,16:,52),istart,icount,mlpt,istat)

         call readvar(filen,'cropout(:,13:,51)',cropout(:,13:,51),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,14:,51)',cropout(:,14:,51),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,15:,51)',cropout(:,15:,51),istart,icount,mlpt,istat)
         call readvar(filen,'cropout(:,16:,51)',cropout(:,16:,51),istart,icount,mlpt,istat)
       endif

! only in the case of
!  - dynamic vegetation restart
!  -  fixed vegetation using the ccm3 climate
! In the other cases, existence arrays are based on climatology
! (inland_main_offline calls climanl that calls existence)
      if ((isimveg .gt. 0) .or. (ccmexist .gt. 0)) then
         do i = lbeg, lend
             tcmin(i) = tc(i) + deltat(i)
          end do
         call existence
      endif

!     return without error
!     --------------------
      return

!     Got error - exit
!     ----------------
9999  write(*,*) 'ERROR in restart, '//varname
      stop 1

end subroutine restartib
