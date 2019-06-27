#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine soilbgc (spin,kpti,kptj)
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control, only: iyear, iyear0, imonth, iday, nspinsoil, &
                                spinmax, ipointout,jday
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comatm
      use inland_comcrop
      use inland_comnitr
      use inland_combgc
      use inland_compft
      use inland_comtex

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer spin,      &   ! # of times soilbgc has been called in the current day (supplied)
              kpti,      &   ! index of 1st point of little vector in big lpt vector
              kptj           ! index of last point of little vector

! local variables
      integer i, j, k
!
      real*8 totts,    &      ! 1/ndaypy
             fracll,   &      ! lignin fraction of leaves 
             fracls,   &      ! structural fraction of leaves 
             fraclm,   &      ! metabolic fraction of leaves 
             fracrl,   &      ! lignin fraction of roots
             fracrs,   &      ! structural fraction of roots 
             fracrm,   &      ! metabolic fraction of roots
             fracwl,   &      ! lignin fraction of wood
             fracws,   &      ! structural fraction of wood
!             cnfroot,  &      ! c/n ratio of fine roots - same value as in biogeochem.f
             fracwm           ! metabolic fraction of wood 
       
      real*8 outclm(lbeg:lend), &   ! c leaving leaf metabolic pool 
             outcls(lbeg:lend), &   ! c leaving leaf structural pool 
             outcll(lbeg:lend), &   ! c leaving leaf lignin pool
             outcrm(lbeg:lend), &   ! c leaving root metabolic pool 
             outcrs(lbeg:lend), &   ! c leaving root structural pool 
             outcrl(lbeg:lend), &   ! c leaving root lignin pool
             outcwm(lbeg:lend), &   ! c leaving woody metabolic carbon pool
             outcws(lbeg:lend), &   ! c leaving woody structural carbon pool
             outcwl(lbeg:lend), &   ! c leaving woody lignin carbon pool
             outcsb(lbeg:lend), &   ! flow of passive c to biomass
             outcps(lbeg:lend), &   ! flow of protected om to passive pool 
             outcns(lbeg:lend), &   ! flow of non-protected om to passive pool
             outcnb(lbeg:lend), &   ! flow of non-protected om to biomass 
             outcpb(lbeg:lend), &   ! flow of protected om to biomass
             outcbp(lbeg:lend), &   ! c leaving protected biomass pool  
             outcbn(lbeg:lend), &   ! c leaving non-protected biomass pool
             totc(lbeg:lend)        ! total c in soil

      real*8 dbdt(lbeg:lend),    &  ! change of c in biomass pools with time 
             dcndt(lbeg:lend),   &  ! change of c in non-protected om with time
             dcpdt(lbeg:lend),   &  ! change of c in protected om with time
             dcsdt(lbeg:lend),   &  ! change of c in passive om with time
             netmin(lbeg:lend),  &  ! net nitrogen mineralization
             nbiors(lbeg:lend),nbiols(lbeg:lend),nbiows(lbeg:lend),      & 
             nbiowm(lbeg:lend),nbiolm(lbeg:lend),nbiorm(lbeg:lend),      & 
             nbioslon(lbeg:lend),nbioslop(lbeg:lend),nbiopas(lbeg:lend), & 
             nminrs(lbeg:lend),nminls(lbeg:lend),nminws(lbeg:lend),      &
             nminwm(lbeg:lend),nminlm(lbeg:lend),nminrm(lbeg:lend),      &
             nminslon(lbeg:lend),nminslop(lbeg:lend),nminpas(lbeg:lend), &
             nrelps(lbeg:lend),nrelns(lbeg:lend),nrelbn(lbeg:lend),      &
             nrelbp(lbeg:lend),nrelll(lbeg:lend),nrelrl(lbeg:lend),      &
             nrelwl(lbeg:lend),ymintot(lbeg:lend),yminmic(lbeg:lend)

! nitrogen in litter and soil pools
      real*8 nlitlm(lbeg:lend),nlitls(lbeg:lend),nlitll(lbeg:lend),      &
             nlitrm(lbeg:lend),nlitrs(lbeg:lend),nlitrl(lbeg:lend),      &
             nlitwm(lbeg:lend),nlitws(lbeg:lend),nlitwl(lbeg:lend),      &
             nsoislop(lbeg:lend),nsoipas(lbeg:lend),nsoislon(lbeg:lend)

! variables controlling constraints on microbial biomass 
      real*8 cmicn(lbeg:lend),cmicp(lbeg:lend),cmicmx(lbeg:lend)

! variables controlling leaching, calculating co2 respiration and n deposition
      real*8 totcbegin(lbeg:lend),totcend(lbeg:lend),totcin(lbeg:lend)

      real*8 fleach

! variables added to do daily time series of some values
      integer kk, textcls, lmin
!
! variables dealing with soil texture and algorithms
      integer msand,mclay
!
      real*8 fsand,fclay,carfrac,texfact,fbpom,rdepth

! total timesteps (daily) used to divide litterfall into daily fractions 
      totts=1./float(ndaypy)

!-----------------------------------------------------------------------
! specific maximum decay rate or growth constants; rates are per day
! constants are taken from Parton et al., 1987 and Verberne et al., 1990
! and special issue of Geoderma (comparison of 9 organic matter models) in Dec. 1997
!
! leaching parameterization was changed to agree with field data,
! this caused a changing of the below constants.  
!
! approximate factors for Verberne et al. model where efficiencies are 100%
! for some of the transformations: one problem was that their rate constants were
! based on 25C, and our modifying functions are based on 15 C...thus the rate constants
! are somewhat smaller compared to the Verberne et al. (1990) model parameters
! rates are based on a daily decomposition timestep (per day)
!-----------------------------------------------------------------------

! leaf litter constants
      klm = klm*kfactor       !dpm leaf --> microbial biomass
      kls = kls*kfactor       !spm leaf --> microbial biomass
      kll = kll*kfactor       !rpm leaf --> non or protected om

! root litter constants
      krm = krm*kfactor       !dpm root --> microbial biomass
      krs = krs*kfactor       !spm root --> microbial biomass
      krl = krl*kfactor       !rpm root --> non or protected om 

! woody litter constants
      kwm = kwm*kfactor       !dpm wood --> microbial biomass
      kws = kws*kfactor       !spm wood --> microbial biomass
      kwl = kwl*kfactor       !rpm wood --> non or protected om 

! biomass constants
      kbn = kbn*kfactor       !biomass --> non protected organic matter 
      kbp = kbp*kfactor       !biomass --> protected organic matter

! slow and passive c pools
      knb = knb*kfactor       !non protected om --> biomass
      kns = kns*kfactor       !non protected om --> stablized om
      kpb = kpb*kfactor       !protected om     --> biomass
      kps = kps*kfactor       !protected om     --> stablized om
      ksb = ksb*kfactor       !stablized om     --> biomass

! ---------------------------------------------------------------------
!  yield (efficiency) with which microbes gain biomass from c source
!  the rest is driven off as co2 respiration (microbial respiration)
!  all of the respiration produced by microbes is assumed to leave
!  the soil profile over the course of a year
!  taken primarily from the models of Verberne and CENTURY
! ---------------------------------------------------------------------
!      ylm = 0.4       ! metabolic material efficiencies
!      yrm = 0.4
!      ywm = 0.4
!      yls = 0.3       ! structural efficiencies
!      yrs = 0.3
!      yws = 0.3

!      yll = 1.0       ! resistant fraction
!      yrl = 1.0 
!      ywl = 1.0 
!      ybn = 1.0       ! biomass       --> non-protected pool
!      ybp = 1.0       ! biomass       --> protected pool
!      yps = 1.0       ! protected     --> passive
!      yns = 1.0       ! non-protected --> passive

!      ysb = 0.20      ! passive pool  --> biomass
!      ypb = 0.20      ! protected     --> biomass
!      ynb = 0.25      ! non-protected --> biomass

! -------------------------------------------------------------------
! split of lignified litter material between protected/non-protected
! slow OM pools
! -------------------------------------------------------------------
!      lig_frac = 0.50 
!
! -------------------------------------------------------------------
! protected biomass as a fraction of total soil organic carbon
! from Verberne et al., 1990
! -------------------------------------------------------------------
!      fbsom = 0.017
!
! ---------------------------------------------------------------------
! (effac) --> efficiency of microbial biomass reincorporated
! into biomass pool.(from NCSOIL parameterizations; Molina et al., 1983)
! ---------------------------------------------------------------------
!      effac = 0.40 
!
! ---------------------------------------------------------------------
! define C:N ratios of substrate pools and biomass
! metabolic, structural, and lignin are for Leaves and roots
! values from Parton et al., 1987 and Whitmore and Parry, 1988
! index: 1 - biomass, 2 - passive pool, 3- slow protected c,
! 4 - slow carbon, non-protected, 5 - resistant, 6 - structural plant
! leaf and root litter, 7 - metabolic plant and root litter, 
! 8- woody biomass
! ---------------------------------------------------------------------
!      cnr(1)  = 8.0       !c:n ratio of microbial biomass
!      cnr(2)  = 15.0      !c:n ratio of passive soil carbon
!      cnr(3)  = 10.0      !c:n ratio of protected slow soil carbon
!      cnr(4)  = 15.0      !c:n ratio of non-protected slow soil C
!      cnr(5)  = 100.0     !c:n ratio of resistant litter lignin
!      cnr(6)  = 150.0     !c:n ratio of structural plant litter
!      cnr(7)  = 6.0       !c:n ratio of metabolic plant litter
!      cnr(8)  = 250.0     !c:n Ratio of woody components
!

! ---------------------------------------------------------------------
! calculate the fraction of wood, roots and leaves that are structural,
! decomposable, and resistant based on equations presented in Verberne
! model discussion (Geoderma, December 1997 special issue).  fmax is the
! maximum fraction allowed in resistant fraction, rconst is a constant
! defined as 1200.  The cnratio of each plant part has to be less than
! the value of structural defined above (i.e. 150) otherwise the equations
! are unstable...thus the wood litter pool value for cnr(6) is substituted
! with a value higher than that for cnwood (i.e. 250).  this is 
! insignificant for wood since 97% is structural anyways.
!
! ** NOTE ******** 
! Would like to incorporate different C:N ratios of residue/roots for
! different biome types based on literature search
! average c:n ratio would be based on litter inputs from each pft
! ****************
! ---------------------------------------------------------------------
!
! equations were changed on 1-26-99 for erratum in literature (Whitmore
! et al. 1997) which had an error in equations to split litterfall into
! the correct three fractions
!      fmax   = 0.45
!      rconst = 1200.0
!      cnwood = 200.0      ! average c:n ratio for woody debris

! leaf litter 
      fracll = fmax * (cnleaf**2)/(rconst + cnleaf**2)
      fracls = (1./cnleaf - fracll/cnr(5) - (1.-fracll)/cnr(7))/     &
               (1./cnr(6) - 1./cnr(7))

      if(isimagro .eq. 0) then
         fraclm = 1.0 - fracll - fracls
      else
         fraclm = max(0.0, 1.0 - fracll - fracls)
      endif

! root litter
         fracrl = fmax * (cnroot**2)/(rconst + cnroot**2)
         fracrs = (1./cnroot - fracrl/cnr(5) - (1.-fracrl)/cnr(7))/&
                  (1./cnr(6) - 1./cnr(7))

         fracrm = 1.0 - fracrl - fracrs

! wood litter
      fracwl = fmax * (cnwood**2)/(rconst + cnwood**2)
      fracws = (1./cnwood - fracwl/cnr(5) - (1.-fracwl)/cnr(7))/     &
               (1./cnr(8) - 1./cnr(7))
      fracwm = 1.0 - fracwl - fracws

      do 100 i = kpti, kptj 

     if(isimagro .gt. 0) then

! check to see if crops are planted - if so, modify c:n ratios 
! crop residue is input on the harvest date


          do 80 j = scpft, ecpft
           if (exist(i,j) .eq. 1 .and. cntops(i,j) .ne. 0 .and. &
               harvdate(i,j) .eq. jday) then
                  cnleaf  = max(40.0, cntops(i,j))   ! calculated in crop residue at harvest
                  cnroot = max(60.0, cnrootvec(i,j))   ! calculated in crop residue at harvest 
                  totts   = 1.            ! pulse input for crops residue
	   elseif(exist(i,16).eq.1.and.croplive(i,16).eq.1) then
                  cnleaf  = max(40.0, cntops(i,j))   ! calculated in crop residue at harvest
                  cnroot = max(60.0, cnrootvec(i,j))   ! calculated in crop residue at harvest 
                  totts   = 1. !put fall daily.
           endif

 80       continue

! leaf litter 
      fracll = fmax * (cnleaf**2)/(rconst + cnleaf**2)
      fracls = (1./cnleaf - fracll/cnr(5) - (1.-fracll)/cnr(7))/     &
               (1./cnr(6) - 1./cnr(7))
      if(isimagro .eq. 0) then
         fraclm = 1.0 - fracll - fracls
      else
         fraclm = max(0.0, 1.0 - fracll - fracls)
      endif
! root litter
         fracrl = fmax * (cnroot**2)/(rconst + cnroot**2)
         fracrs = (1./cnroot - fracrl/cnr(5) - (1.-fracrl)/cnr(7))/&
                  (1./cnr(6) - 1./cnr(7))
      fracrm = 1.0 - fracrl - fracrs

! wood litter
      fracwl = fmax * (cnwood**2)/(rconst + cnwood**2)
      fracws = (1./cnwood - fracwl/cnr(5) - (1.-fracwl)/cnr(7))/     &
               (1./cnr(8) - 1./cnr(7))
      fracwm = 1.0 - fracwl - fracws

     endif ! check for crop existence


! ---------------------------------------------------------------------
! fraction of decomposing microbial biomass into protected organic
! matter; taken from the model of Verberne et al., 1990
! this is the proportion of decomposing dead microbial biomass that
! is transferred to a protected pool vs. a non-protected pool
! related to the clay content of the soil. in sandy soils, fbpom = 0.3,
! whereas in clay soils fbpom = 0.7.  created a linear function based
! on clay fraction of soil to adjust according to amount of clay in
! the top 1 m of profile (weighted average according to depth of each
! layer)
!
! also take care of calculation of texfact, which is a leaching
! parameter based on the average sand fraction of the top 1 m of
! soil
! ---------------------------------------------------------------------
       rdepth = 0.0
       do kk = 1, nslaym
        rdepth     = rdepth + hsoi(kk)
       enddo 
         
         rdepth   = 1./rdepth
         carfrac    = 0.0
         texfact  = 0.0 

        do k = 1, nslaym
            !Fontes: infilensoilayer is equal 0 on single_point
          if (k.le.infilensoilayer.or.infilensoilayer.eq.0) then
            msand = nint(sand(i,k))
            mclay = nint(clay(i,k))
          else
            msand = nint(sand(i,infilensoilayer))
            mclay = nint(clay(i,infilensoilayer))
          endif
	enddo

        do 90 kk = 1, nslaym                    ! top 1 m of soil -- 8 layers
            msand    = nint(sand(i,kk)) 
            mclay    = nint(clay(i,kk)) 
            fclay    = 0.01 * mclay
            fsand    = 0.01 * msand  
            carfrac    = carfrac   + fclay * hsoi(kk)
            texfact  = texfact + fsand * hsoi(kk)
90      continue

         carfrac   = carfrac   * rdepth
         texfact = texfact * rdepth

! if carfrac is greater than 0.4, set fbpom = 0.7, if carfrac is less
! than 0.17, set fbpom = 0.30 (sandy soil)
!
!       fbpom = min(max(0.3, carfrac/0.4 * 0.7),0.7)      
         fbpom = 0.50

! ------------------------------------------------------------------------
! total soil carbon initialized to 0 at beginning of model run
! used in calculation of soil co2 respiration from microbial decomposition 
! ------------------------------------------------------------------------
         if (iday .eq. 1 .and. imonth .eq. 1 .and. iyear .eq. iyear0) then
            totcbegin(i) = 0.0
            storedn(i)   = 0.0 
         endif

! ------------------------------------------------------------------------
! initialize yearly summation of net mineralization and co2 respiration
! to 0 at beginning of each year; because these quantities are usually 
! reported on a yearly basis, we wish to do the same in the model so we
! can compare easily with the data.
! ------------------------------------------------------------------------
         if (iday .eq. 1 .and. imonth .eq. 1) then
            yrleach(i) = 0.0
            cleach(i)  = 0.0
            ynleach(i) = 0.0
            ymintot(i) = 0.0
            yminmic(i) = 0.0
         endif
        
! determine amount of substrate available to microbial growth
!
! calculate the total amount of litterfall entering soil(C)
         totcin(i) = falll(i)*totts + fallr(i)*totts + fallw(i)*totts

! calculate the current total amount of carbon at each grid cell
         totc(i) = clitlm(i) + clitls(i) + clitrm(i) + clitrs(i) +      &
                   clitwm(i) + clitws(i) + csoislop(i) + csoislon(i) +    &
                   csoipas(i) + totcmic(i) + clitll(i) + clitrl(i) + clitwl(i)

! beginning amount of soil C at each timestep (used for respiration
! calculation)
         totcbegin(i) = totc(i)

! ------------------------------------------------------------------------
! split current amount of total soil microbes
! maximum amount of biomass is a function of the total soil C
! from Verberne et al., 1990
! ------------------------------------------------------------------------
!       totcmic(i) = cmicp(i) + cmicn(i)
         cmicmx(i) = fbsom * totc(i) 

! calculate the amount of protected and unprotected biomass
         if (totcmic(i) .ge. cmicmx(i)) then
            cmicp(i) = cmicmx(i)
            cmicn(i) = totcmic(i) - cmicmx(i)
         else
            cmicn(i) = 0.0
            cmicp(i) = totcmic(i)
         endif

! ---------------------------------------------------------------
! litter pools 
!
! add in the amount of litterfall, and root turnover
! ---------------------------------------------------------------
         clitlm(i) = clitlm(i) + (fraclm * falll(i)*totts)  
         clitls(i) = clitls(i) + (fracls * falll(i)*totts)  
         clitll(i) = clitll(i) + (fracll * falll(i)*totts)  
         clitrm(i) = clitrm(i) + (fracrm * fallr(i)*totts)  
         clitrs(i) = clitrs(i) + (fracrs * fallr(i)*totts)  
         clitrl(i) = clitrl(i) + (fracrl * fallr(i)*totts)  
         clitwm(i) = clitwm(i) + (fracwm * fallw(i)*totts)  
         clitws(i) = clitws(i) + (fracws * fallw(i)*totts)  
         clitwl(i) = clitwl(i) + (fracwl * fallw(i)*totts) 

! ---------------------------------------------------------------
! calculate microbial growth rates based on available C sources
! to microbes (substrate : litter, C in slow, passive pools)
! the amount of biomass added cannot be larger than the amount of
! available carbon from substrates and other pools at this point.
! ---------------------------------------------------------------
         outcrs(i) = min(decomps(i) * krs * clitrs(i),clitrs(i))
         outcws(i) = min(decompl(i) * kws * clitws(i),clitws(i))
         outcls(i) = min(decompl(i) * kls * clitls(i),clitls(i))
         outclm(i) = min(decompl(i) * klm * clitlm(i),clitlm(i))
         outcrm(i) = min(decomps(i) * krm * clitrm(i),clitrm(i))
         outcwm(i) = min(decompl(i) * kwm * clitwm(i),clitwm(i))
         outcnb(i) = min(decomps(i) * knb * csoislon(i),csoislon(i))
         outcpb(i) = min(decomps(i) * kpb * csoislop(i),csoislop(i))
         outcsb(i) = min(decomps(i) * ksb * csoipas(i),csoipas(i))

! ---------------------------------------------------------------
! calculate turnover of microbial biomass
! two disctinct pools: one with rapid turnover, and one with slow
! turnover rate
! ---------------------------------------------------------------
         outcbp(i) = min(kbp * cmicp(i),cmicp(i))
         outcbn(i) = min(kbn * cmicn(i),cmicn(i))

! ---------------------------------------------------------------------
! recycle microbes back to respective microbial pools based on effac as
! discussed in NCSOIL model from Molina et al., 1983
! ---------------------------------------------------------------------
         outcbp(i) = outcbp(i) *  effac
         outcbn(i) = outcbn(i) *  effac

! -------------------------------------------------------------------------
! have to adjust inputs into microbial pool for the slow
! and passive carbon amounts that are leaving their respective
! pools at an increased rate during the spinup procedure.
! these values should be decreased by the respective spinup factors
! because the microbial pools will otherwise become larger without
! scientific reason due to the spinup relationships used.
! 3 main pools: outcpb, outcnb, outcsb
! -------------------------------------------------------------------------
         dbdt(i) = outcrs(i) * yrs + outcws(i) * yws +     &
                   outcls(i) * yls + outclm(i) * ylm +     &
                   outcrm(i) * yrm + outcwm(i) * ywm +     &
                   outcnb(i) * ynb + outcpb(i) * ypb +     &
                   outcsb(i) * ysb - outcbp(i) - outcbn(i)

! -------------------------------------------------------------------------
! change in non-protected organic matter from growth in microbial
! biomass, lignin input, and stablized organic matter pool
! the flow out of the pool from its decomposition is always less
! the yield--which is factored into the pool it is flowing into
! -------------------------------------------------------------------------
#ifdef SINGLE_POINT_MODEL /* This problem probably only occurs on 0D version! */
! FIXME: several variables can become too small and make number converge to zero
!        For example, for a KM83 simulation, starting in 1-jan-1992, the follow-
!        ing variables become too small (1.0E-30) in:
!        outclm: 8-sep-1992: c leaving leaf metabolic pool
!        outcrm: 11-dec-1992: c leaving root metabolic pool
!        outcwm: 27-jan-1993: c leaving wood metabolic pool
!        outcll: 12-mar-1999: c leaving leaf lignin pool (in some posterior
!                steps this oscillates between more and less than 1E-30, but
!                converges to smaller than 1E-30 after somw while)
!        outcls: 06-sep-2001: c leaving leaf structural pool
#endif /* SINGLE_POINT_MODEL */
         outcll(i) = min(decompl(i) * kll * clitll(i),clitll(i))
         outcrl(i) = min(decomps(i) * krl * clitrl(i),clitrl(i))
         outcwl(i) = min(decompl(i) * kwl * clitwl(i),clitwl(i))
         outcns(i) = min(decomps(i) * kns * csoislon(i),csoislon(i))

! ------------------------------------------------------------ 
! the lig_fac factor only applies to lignin content...half goes to
! protected slow OM, and half goes to non protected slow OM
! ------------------------------------------------------------
         dcndt(i) = (lig_frac * (outcll(i) * yll + outcrl(i) * yrl +     &
                    outcwl(i) * ywl) + (1. - fbpom) * (ybn * outcbn(i) + &
                    ybp * outcbp(i))) - outcnb(i) - outcns(i)

! ------------------------------------------------------------
! change in protected organic matter from growth in microbial 
! biomass, lignin input, and stablized organic matter pool
! ------------------------------------------------------------
         outcps(i) = min(decomps(i) * kps * csoislop(i), csoislop(i))

! ------------------------------------------------------------
! the lig_frac factor only applies to lignin content...half goes to
! protected slow OM, and half goes to non protected slow OM
! ------------------------------------------------------------
         dcpdt(i) = (lig_frac * (outcll(i)*yll+outcrl(i)*yrl +    &
                    outcwl(i) * ywl) + fbpom * (ybn * outcbn(i) + &
                    ybp * outcbp(i))) - outcpb(i) - outcps(i)

! ----------------------------------------------------------------------
! change in stablized organic matter (passive pool) from growth
! in microbial biomass, and changes in protected and unprotected
! SOM
!
! add a loss of C due to leaching out of the profile, based
! on approximation of CENTURY model below 1 m in depth
! based on water in the profile, and texture of soil
! tuned to known outputs or leaching that has been measured in the field
! at Arlington-WI (Courtesy K. Brye, MS) and applied to the global scale
! on average, this calibration yields about 10-50 Kg C ha-1 yr-1 leaching
! depending on C in soil...will need to be tied to an amount of water
! flowing through the profile based upon precipitation eventually
! ----------------------------------------------------------------------
!         h20    = 0.30e-03

! h20 is a constant relating to the flow of water through the top 1 m of the
! profile 
! use texfact -- the % sand -- or texture factor effect on leaching (see Parton
! et al. (1991) calculated from the average sand content of top 1 m of soil
! in the model
         fleach = h20/18.0 * (0.01 + 0.04 * texfact)

! --------------------------------------------------------------------
! change in passive organic carbon pool
! ---------------------------------------------------------------------
         dcsdt(i) = ((yns * outcns(i)) + (yps * outcps(i))) -      &
                     outcsb(i) -  (fleach * csoipas(i))
         cleach(i) = fleach * csoipas(i) + fleach * csoislop(i) +  &
                     fleach * csoislon(i)
         ynleach(i) = ynleach(i) + fleach * csoipas(i)/cnr(2) +  &
                      fleach * csoislop(i)/cnr(3) +              &
                      fleach * csoislon(i)/cnr(4)

! update slow pools of carbon for leaching losses
         dcndt(i) = dcndt(i) - fleach * csoislon(i)              
         dcpdt(i) = dcpdt(i) - fleach * csoislop(i)
         if (spin .eq. spinmax) then
            yrleach(i) =  cleach(i) + yrleach(i)
         endif

! ---------------------------------------------------------------------
! calculate the amount of net N mineralization or immobilization
! ---------------------------------------------------------------------

! uptake of n by growth of microbial biomass
!
! immobilized n used for requirements of microbial growth 
! is based on flow of carbon and the difference of C/N ratio of
! the microbes and their efficiency versus the C/N ratio of the
! material that is being decomposed 
!
! ------------------------------
! structural root decomposition 
! ------------------------------
         if (yrs/cnr(1) .gt. 1./cnr(6)) then
            nbiors(i) = (1./cnr(6) - yrs/cnr(1))* outcrs(i)
            nminrs(i) = 0.0
         else
            nminrs(i) = (1./cnr(6) - yrs/cnr(1))* outcrs(i)
            nbiors(i) = 0.0
         endif

! ------------------------------
! structural leaf decomposition
! ------------------------------
         if (yls/cnr(1) .gt. 1./cnr(6)) then
            nbiols(i) = (1./cnr(6) - yls/cnr(1))* outcls(i)
            nminls(i) = 0.0
         else
            nminls(i) = (1./cnr(6) - yls/cnr(1))* outcls(i)
            nbiols(i) = 0.0
         endif

! ------------------------------
! structural wood decomposition
! ------------------------------
         if (yws/cnr(1) .gt. 1./cnr(8)) then
            nbiows(i) = (1./cnr(8) - yws/cnr(1))* outcws(i)
            nminws(i) = 0.0
         else
            nminws(i) = (1./cnr(8) - yws/cnr(1))*outcws(i)
            nbiows(i) = 0.0
         endif

! ------------------------------
! metabolic wood decomposition
! ------------------------------
         if (ywm/cnr(1) .gt. 1./cnr(8)) then
            nbiowm(i) = (1./cnr(8) - ywm/cnr(1))* outcwm(i)
            nminwm(i) = 0.0
         else
            nminwm(i) = (1./cnr(8) - ywm/cnr(1))* outcwm(i)
            nbiowm(i) = 0.0
         endif

! ------------------------------
! metabolic leaf decomposition
! ------------------------------
         if (ylm/cnr(1) .gt. 1./cnr(7)) then
           nbiolm(i) = (1./cnr(7) - ylm/cnr(1))* outclm(i)
           nminlm(i) = 0.0
         else
            nminlm(i) = (1./cnr(7) - ylm/cnr(1))* outclm(i)
            nbiolm(i) = 0.0
         endif

! ------------------------------
! metabolic root decomposition
! ------------------------------
         if (yrm/cnr(1) .gt. 1./cnr(7)) then
            nbiorm(i) = (1./cnr(7) - yrm/cnr(1))* outcrm(i)
            nminrm(i) = 0.0
         else
            nminrm(i) = (1./cnr(7) - yrm/cnr(1))* outcrm(i)
            nbiorm(i) = 0.0
         endif

! ----------------------------------------------
! non-protected organic matter decomposition
! ----------------------------------------------
         if (ynb/cnr(1) .gt. 1./cnr(4)) then
            nbioslon(i) = (1./cnr(4) - ynb/cnr(1))* outcnb(i)
            nminslon(i) = 0.0
         else
            nminslon(i) = (1./cnr(4) - ynb/cnr(1))* outcnb(i)
            nbioslon(i) = 0.0
         endif

! ----------------------------------------------
! protected organic matter decomposition
! ----------------------------------------------
         if (ypb/cnr(1) .gt. 1./cnr(3)) then
            nbioslop(i) = (1./cnr(3) - ypb/cnr(1))* outcpb(i)
            nminslop(i) = 0.0
         else
            nminslop(i) = (1./cnr(3) - ypb/cnr(1))* outcpb(i)
            nbioslop(i) = 0.0
         endif

! ----------------------------------------------
! stablized organic matter decomposition
! ----------------------------------------------
         if (ysb/cnr(1) .gt. 1./cnr(2)) then
            nbiopas(i) = (1./cnr(2) - ysb/cnr(1))* outcsb(i)
            nminpas(i) = 0.0
         else
            nminpas(i) = (1./cnr(2) - ysb/cnr(1))* outcsb(i)
            nbiopas(i) = 0.0
         endif

! ----------------------------------------------
! total immobilized N used for biomass growth
! ----------------------------------------------
         totimm(i) = nbiors(i) + nbiols(i) + nbiows(i) + nbiowm(i)     &
                   + nbiolm(i) + nbiorm(i) + nbioslon(i) + nbioslop(i) &
                   + nbiopas(i)

! -----------------------------------------------------------------------------
! gross amount of N mineralized by decomposition of C by microbial biomass
! assume that N is attached to the flow of C by the C/N ratio of the substrate
! also assume that the amount of N attached to CO2 that is respired is also
! mineralized (i.e. the amount of N mineralized is related to the total outflow
! of carbon, and not the efficiency or yield)..see Parton et al., 1987
! -----------------------------------------------------------------------------
         totmin(i) = nminrs(i) + nminls(i) + nminws(i) + nminwm(i)      &
                   + nminlm(i) + nminrm(i) + nminslon(i) + nminslop(i)  &
                   + nminpas(i)

! -----------------------------------------------------------------------------
! when carbon is transferred from one pool to another, each pool has a distinct
! C:N ratio.  In the case of pools where carbon is moving from the pool to 
! the microbial biomass (used for growth/assimilation), net mineralization
! takes place (N is released) after the requirements of building the biomass
! are met.  In the cases of other transformations of C, N is not conserved
! if it follows from one pool to another which has a different C:N ratio;
! either N is released or is needed to make the transformation and keep N
! conserved in the model. 
!
! other calculations of either N release or immobilization to keep track of
! the budget
         nrelps(i) = outcps(i) * (1./cnr(3) - 1./cnr(2))
         nrelns(i) = outcns(i) * (1./cnr(4) - 1./cnr(2))
         nrelbn(i) = (1.-fbpom) * outcbn(i) * (1./cnr(1) - 1./cnr(4)) +  &
                     (1.-fbpom) * outcbp(i) * (1./cnr(1) - 1./cnr(4))
         nrelbp(i) = fbpom * outcbp(i) * (1./cnr(1) - 1./cnr(3)) +       &
                     fbpom * outcbn(i) * (1./cnr(1) - 1./cnr(3))
         nrelll(i) = lig_frac * outcll(i) * (1./cnr(5) - 1./cnr(3)) +    &
                     lig_frac * outcll(i) * (1./cnr(5) - 1./cnr(4))
         nrelrl(i) = lig_frac * outcrl(i) * (1./cnr(5) - 1./cnr(3)) +    &
                     lig_frac * outcrl(i) * (1./cnr(5) - 1./cnr(4))
         nrelwl(i) = lig_frac * outcwl(i) * (1./cnr(5) - 1./cnr(3)) +    &
                     lig_frac * outcwl(i) * (1./cnr(5) - 1./cnr(4))

         totnrel(i) = nrelps(i) + nrelns(i) + nrelbn(i) +                &
                      nrelbp(i) + nrelll(i) + nrelrl(i) + nrelwl(i)

! -----------------------------------------------------------------------------
! calculate whether net mineralization or immobilization occurs
! on a grid cell basis -- tnmin is an instantaneous value for each time step
! it is passed along to stats to calculate, daily, monthly and annual totals
! of nitrogen mineralization
! this is for mineralization/immobilization that is directly related to 
! microbial processes (oxidation of carbon)
!
! the value of totnrel(i) would need to be added to complete the budget
! of N in the model. Because it can add/subtract a certain amount of N
! from the amount of net mineralization.  However, these transformations
! are not directly related to microbial decomposition, so do we add them
! into the value or not?
! -----------------------------------------------------------------------------
         netmin(i) = totmin(i) + totimm(i) + totnrel(i) 
         if (netmin(i) .gt. 0.0) then
            tnmin(i) = netmin(i)
         else
            tnmin(i) = 0.0 
         endif

! convert value of tnmin of Kg-N/m2/dtime to mole-N/s
! based on N = .014 Kg/mole -- divide by the number of seconds in daily timestep
!
            tnmin(i)   = tnmin(i)/(86400.   * 0.014)
          if(isimagro .gt. 0)then
            totmin(i)  = totmin(i)/(86400.  * 0.014)
            totimm(i)  = totimm(i)/(86400.  * 0.014)
            totnrel(i) = totnrel(i)/(86400. * 0.014)
          endif
!
! ---------------------------------------------------
! update soil c pools for transformations of c and n
! ---------------------------------------------------
         totcmic(i)  = max(totcmic(i)  + dbdt(i), dble(0.0))
         csoislon(i) = max(csoislon(i) + dcndt(i),dble(0.0))
         csoislop(i) = max(csoislop(i) + dcpdt(i),dble(0.0))
         csoipas(i)  = max(csoipas(i)  + dcsdt(i),dble(0.0))
         clitlm(i)   = max(clitlm(i)  - outclm(i),dble(0.0))
         clitls(i)   = max(clitls(i)  - outcls(i),dble(0.0))
         clitll(i)   = max(clitll(i)  - outcll(i),dble(0.0))
         clitrm(i)   = max(clitrm(i)  - outcrm(i),dble(0.0))
         clitrs(i)   = max(clitrs(i)  - outcrs(i),dble(0.0))
         clitrl(i)   = max(clitrl(i)  - outcrl(i),dble(0.0))
         clitwm(i)   = max(clitwm(i)  - outcwm(i),dble(0.0))
         clitws(i)   = max(clitws(i)  - outcws(i),dble(0.0))
         clitwl(i)   = max(clitwl(i)  - outcwl(i),dble(0.0))

! -----------------------------------------------------------
! update soil n pools based on c:n ratios of each pool
! this approach is assuming that the c:n ratios are remaining
! constant through the simulation. flow of nitrogen is attached
! to carbon 
! -----------------------------------------------------------
         totnmic(i)  = totcmic(i) /cnr(1)
         nsoislon(i) = csoislon(i)/cnr(4)
         nsoislop(i) = csoislop(i)/cnr(3)
         nsoipas(i)  = csoipas(i) /cnr(2)
         nlitlm(i)   = clitlm(i)  /cnr(7)
         nlitls(i)   = clitls(i)  /cnr(6)
         nlitll(i)   = clitll(i)  /cnr(5)
         nlitrm(i)   = clitrm(i)  /cnr(7)
         nlitrs(i)   = clitrs(i)  /cnr(6)
         nlitrl(i)   = clitrl(i)  /cnr(5)
         nlitwm(i)   = clitwm(i)  /cnr(8)
         nlitws(i)   = clitws(i)  /cnr(8)
         nlitwl(i)   = clitwl(i)  /cnr(8)

! total above and belowground litter
         totlit(i) = clitlm(i) + clitls(i) + clitll(i) +  &
                     clitrm(i) + clitrs(i) + clitrl(i) +  &
                     clitwm(i) + clitws(i) + clitwl(i)

! sum total aboveground litter (leaves and wood)
         totalit(i) = clitlm(i) + clitls(i) + clitwm(i) +  &
                      clitll(i) + clitws(i) + clitwl(i)

! sum total belowground litter (roots) 
         totrlit(i) = clitrm(i) + clitrs(i) + clitrl(i)

! determine total soil carbon amounts (densities are to 1 m depth; Kg/m-2)
         totcsoi(i) = csoipas(i) + csoislop(i) + totcmic(i) + csoislon(i)

! calculate total amount of litterfall occurring (total for year)
         totfall(i) = falll(i) + fallr(i) + fallw(i)
      
         if (spin .eq. spinmax) then
!          if (ipointout .eq. 1 .and. myid .eq. 0) then
            if (ipointout .eq. 1) then
               if (i .eq. 1001) then
                  write(122,'(4f10.5)')                &
                  totrlit(i)*1e2, totlit(i), totalit(i), totcsoi(i)
!                call flush(122)
               endif

               if (i .eq. 662) then
                  write(1122,'(4f10.5)')                 &
                  totrlit(i)*1e2, totlit(i), totalit(i), totcsoi(i)
!                call flush(1122)
               endif

               if(i .eq. 651) then
                  write(2122,'(4f10.5)')                 &
                  totrlit(i)*1e2, totlit(i), totalit(i), totcsoi(i)
!                call flush(2122)
               endif
 
               if (i .eq. 1001) write(*,*)'soilbgc ',totrlit(i),   &
                  totlit(i),totalit(i), totcsoi(i)
               endif
            endif

! nitrogen 
!
! total nitrogen in litter pools (above and belowground)
!
         totnlit(i) = nlitlm(i) + nlitls(i) + nlitrm(i) + nlitrs(i) +  &
                      nlitwm(i) + nlitws(i) + nlitll(i) + nlitrl(i) +  &
                      nlitwl(i)
!
! sum total aboveground litter   (leaves and wood)
!
         totanlit(i) = nlitlm(i) + nlitls(i) + nlitwm(i) +    &
                       nlitll(i) + nlitws(i) + nlitwl(i)
!
! sum total belowground litter  (roots)
!
         totrnlit(i) = nlitrm(i) + nlitrs(i) + nlitrl(i)
!
! total soil nitrogen to 1 m depth (kg-N/m**2)
!  
         totnsoi(i) = nsoislop(i) + nsoislon(i) +              &
                      nsoipas(i)  + totnmic(i) + totnlit(i)

! --------------------------------------------------------------------------
! calculate running sum of yearly net mineralization, and nitrogen in pool
! available to plants for uptake--during spin up period, can only count one
! of the cycles for each timestep--otherwise false additions will result
! values of yearly mineralization are in Kg/m-2
! --------------------------------------------------------------------------
         if (spin .eq. spinmax) then
            storedn(i)  = storedn(i) + tnmin(i)
         endif

! calculate total amount of carbon in soil at end of cycle
! this is used to help calculate the amount of carbon that is respired
! by decomposing microbial biomass
         totcend(i) = totlit(i) + totcsoi(i)

! --------------------------------------------------------------------------
! the amount of co2resp(i) is yearly value and is dependent on the amount
! of c input each year, the amount in each pool at beginning of the year,
! and the amount left in the pool at the end of the year
! along with the amount of root respiration contributing to the flux from
! calculations performed in stats.f
! --------------------------------------------------------------------------
         if (spin .eq. spinmax) then 

! --------------------------------------------------------------------------
! only count the last cycle in the spin-up for co2soi
! when the iyear is less than the nspinsoil value...otherwise
! an amount of CO2 respired will be about 10 times the actual
! value because this routine is called articially 10 extra times
! each time step to spin up the soil carbon
!
! add n-deposition due to rainfall once each year, and
! the amount of N fixed through N-fixers.  These equations
! are based on the annual precip input (cm) and are from
! the CENTURY model...Parton et al., 1987.
! The base equations are in units of (g) N m-2 so have to
! divide by 1000 to put in units of Kg.
!
! deposition now calculated every day. 
! 0.0005753 = 0.21/365
! (-0.0004932 = -0.18/365
! - the original constants (0.21, -0.18) are for the entire year 
! --------------------------------------------------------------------------
         if(isimagro .eq.0 )then
            deposn(i) =(0.0005753  + 0.0028 * ((adrain(i)+adsnow(i))*0.1))*1.e-3
            fixsoin(i)=(-0.0004932 + 0.14  * ((adrain(i)+adsnow(i))*0.1))*1.e-3
         else
            deposn(i) =(0.0005753  + 0.0028 * (precip(i)*0.1))*1.e-3
            fixsoin(i)=(-0.0004932 + 0.14  * (precip(i)*0.1))*1.e-3

! 7-06-05 CJK : modify rates for the US based on Simon's dataset that adjusts
! the rates each year (from actual data - 1940-1999)
!
            if (iyear .ge. 1940 .and. iyear .le. 1999) deposn(i) = deposn(i) * &
                                 0.5 !ndepfact(i,iyear+1-1940) 
            if (iyear .gt. 1999) deposn(i) = deposn(i) * &
                                 0.5 !ndepfact(i,60) 
!
! From Landsberg and Gower, 1997
! changes cjk 5-21-01 : research on rates of asymbiotic nitrogen fixation
! which would be characterized more by IBIS than symbiotic, are about 1-2% of
! the symbiotic potential.  So, adjusting fixsoin by that ratio.
! also maximum allowable annual fixsoin is now 30 kg N ha-1 based on
! Landsberg and Gower, 1997
!
! typical asymbiotic n-fixation by natural vegetation ranges from 0.1 -6
! kg N ha-1 in decidous and coniferous (temperate and boreal forests) to
! 20 in tropical forests.
!
            if (yfixsoin(i) .ge. 0.0030) then
              fixsoin(i) = 0.0 
            else 
              fixsoin(i) = fixsoin(i) * 0.02
            endif
!
            ydeposn(i)   = ydeposn(i)  + deposn(i)
            yfixsoin(i)  = yfixsoin(i) + fixsoin(i)
         
         endif ! check for crop existence
!
! --------------------------------------------------------------------------
! add to the daily total of co2 flux leaving the soil from microbial
! respiration -- instantaneous value for each timestep
! since this subroutine gets called daily...instantaneous fluxes
! the fluxes need to be put on a per second basis, which will be dependent
! on the timestep.  Furthermore, because the biogeochem subroutine does
! not get called each timestep...an approximation for a timestep average
! microbial flux and nmineralization rate will be applied
! --------------------------------------------------------------------------
!
! calculate daily co2 flux due to microbial decomposition
            tco2mic(i) = totcbegin(i) + totcin(i) - totcend(i) - cleach(i)
           
! NOTICE: Specific for the single point version, tco2mic is far behind the
!        difference the original ibis 0D version finds. Where 0D IBIS gets a value
!        of about 2.7E-06, the InLand 0D version, where variables are double
!        precision REAL, it is changed down to 1.57E-15!
! TODO: Analyse this effect in a cumulative/long run.

! convert co2 flux from kg C/day  (seconds in a daily timestep) to mol-C/s
! based on .012 Kg C/mol
            tco2mic(i) = tco2mic(i)/(86400. * .012)
         endif
100   continue

! return to main
      return
end subroutine soilbgc
