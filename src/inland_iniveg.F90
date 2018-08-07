#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine iniveg (isimveg, irestart)
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_combcs
      use inland_comcrop
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_compft
      use inland_combgc
      
      implicit none
!-----------------------------Arguments---------------------------------
! Input arguments
      integer isimveg, &        ! 0:static veg 1: dyn veg 2: dyn veg-cold start
              irestart          ! 0: not a restart run 1: restart run

! local variables
      integer npts,    &
              ideci,   &    ! # deciduous plant functional types (pft)
              ievgr,   &    ! # evergreen pft 
              ishrub,  &    ! # shrub pft 
              igrass,  &    ! # herbaceous pft
              icrop,   &    ! # cropland pft 
              ilower,  &    ! possible # pft for lower canopy
              iupper        ! possible # pft for upper canopy

      integer k,p,j,i     ! loop indices
      integer inveg        ! vegetation type

      real*8 plaievgr,  &  ! potential lai of evergreen trees
             plaideci,  &  ! potential lai of deciduous trees
             plaishrub, &  ! potential lai of shrubs
             plaigrass, &  ! potential lai of grasses
             wood,      &  ! total wood biomas in grid cell
             totdepth,  &  ! total soil depth
             frootnorm1,&  ! normalization factor for Jackson rooting profile,low
             frootnorm2    ! normalization factor for Jackson rooting profile, up

      real*8 depth(nsoilay) ! soil layer depth (cm)

! initialize a few climatic variables needed for vegetation
      if (irestart .eq. 0) then
             agddu(:) = 1000.0
             agddl(:) = 1000.0
             cntops(:,:) = 40.0 
             cnrootvec(:,:) = cnroot
             tnplant(:,:) = 0.0 
             harvdate(:,:) = 999
! initialize the moisture stress factors
         stresstu(:) = 1.0
         stresstl(:) = 1.0
         stressl(:,:) = 1.0
         stressu(:,:) = 1.0

      endif
!
     if(isimagro .eq. 0) woodnorm = 7.5
! Kaiyuan Li for root-water uptake
        stre_tu(:) = 1.0
        stre_tl(:) = 1.0


! initialize running-mean air temperature if not restart
      if (irestart .eq. 0) then

! we are already initializing TD using the weather generator if not coupled
!           a10td(:) = 273.16
           a11soiltd(:) = 273.16
           daylength(:) = 0.
           a3tdmin(:) = 0.

#ifndef SINGLE_POINT_MODEL
         a10td(:) = 273.16
#else /* SINGLE_POINT_MODEL */
         a10td(:) = 293.16
#endif /* SINGLE_POINT_MODEL */

! initialize running-mean values of canopy photosynthesis rates 
! if not restart
         a10ancub(:) = 10.0e-06
         ! It is said a10ancuc is not being used, so its commented out - fzm
         !a10ancuc(i) = 10.0e-06
         a10ancls(:) = 10.0e-06
         a10ancl4(:) = 10.0e-06
         a10ancl3(:) = 10.0e-06

! initialize running-mean values of the scaling parameter
         a10scalparamu(:) = 0.5 * 5.
         a10scalparaml(:) = 0.5 * 5.
         a10daylightu(:) = 5.
         a10daylightl(:) = 5.

! initialize litter fall
         falll(:) = 0.
         fallr(:) = 0.
         fallw(:) = 0.

! initialize this year's growing degree days and temperature of the
! warmest and coldest month for existence of pfts
         gdd0this(:) = 0.
         gdd5this(:) = 0.
         tcthis(:) = 100.
         twthis(:) = - 100.        
      endif

      do 100 i = lbeg, lend
         ievgr  = 0
         ideci  = 0
         ishrub = 0
         igrass = 0
         icrop  = 0
         ilower = 0
         iupper = 0

! determine number of evergreen plant functional types
         if (nint(exist(i,1)).eq.1) ievgr = ievgr + 1
         if (nint(exist(i,3)).eq.1) ievgr = ievgr + 1
         if (nint(exist(i,4)).eq.1) ievgr = ievgr + 1
         if (nint(exist(i,6)).eq.1) ievgr = ievgr + 1

! determine number of deciduous plant functional types
         if (nint(exist(i,2)).eq.1) ideci = ideci + 1
         if (nint(exist(i,5)).eq.1) ideci = ideci + 1
         if (nint(exist(i,7)).eq.1) ideci = ideci + 1
         if (nint(exist(i,8)).eq.1) ideci = ideci + 1

! make sure counter is at least 1 (to avoid division by zero)
         ievgr = max (1, ievgr)
         ideci = max (1, ideci)

! determine number of shrub functional types
         if (nint(exist(i,9)).eq.1)  ishrub = ishrub + 1
         if (nint(exist(i,10)).eq.1) ishrub = ishrub + 1

! determine number of herbaceous plant functional types
         if (nint(exist(i,11)).eq.1) igrass = igrass + 1
         if (nint(exist(i,12)).eq.1) igrass = igrass + 1

! make sure counter is at least 1 (to avoid division by zero)
         ishrub = max (1, ishrub)
         igrass = max (1, igrass)

! total number of possible pfts for each canopy
         iupper = ievgr  + ideci
         ilower = ishrub + igrass

! make sure counter is at least 1 (to avoid division by zero)
         iupper = max (1, iupper)
         ilower = max (1, ilower)
     if(isimagro .gt. 0)then
!
! determine number of crop functional types
         if (nint(exist(i,13)).eq.1) then
            icrop = icrop + 1   ! soybean 
            if(irestart .eq. 0) then
               croplive(i,13) = 0.0
            endif !end irestart
         endif
         if (nint(exist(i,14)).eq.1) then
            icrop = icrop + 1   ! maize
            if(irestart .eq. 0) then
               croplive(i,14) = 0.0
            endif!end irestart
         endif
         if (nint(exist(i,15)).eq.1) then
            icrop = icrop + 1   ! wheat
            if(irestart .eq. 0) then
               croplive(i,15) = 0.0
            endif!end irestart
         endif
         if (nint(exist(i,16)).eq.1) then
            icrop = icrop + 1   ! sugarcane
            if(irestart .eq. 0) then
               croplive(i,16) = 0.0
            endif!end irestart
         endif

! make sure counter is at least 1 (to avoid division by zero)
!
        icrop = max (1, icrop)
        ilower = ilower + icrop
!
        iupper = max (1, iupper)
        ilower = max (1, ilower)
!
     endif ! check for crop existence
! cold start of the vegetation
         if (irestart .eq. 0) then

! ************************************************************************
! case (0) assign vegetation characteristics for static vegetation
! ************************************************************************
!
! and
!
! ************************************************************************
! case (1) assign vegetation characteristics for dynamic vegtation
!          that is initialized with fixed vegetation map
! ************************************************************************
            if (isimveg.eq.0 .or. isimveg.eq.1) then

! translate vegetation type (real) to nearest integer
               inveg = nint (xinveg(i))

! for initialization purposes, set the predicted vegetation type
! to the initial vegetation type
               vegtype0(i) = xinveg(i)

! ---------------------------------------------------
!  1: tropical evergreen forest / woodland
!  2: tropical deciduous forest / woodland
!  3: temperate evergreen broadleaf forest / woodland
!  4: temperate evergreen conifer forest / woodland
!  5: temperate deciduous forest / woodland
!  6: boreal evergreen forest / woodland
!  7: boreal deciduous forest / woodland
!  8: mixed forest / woodland
!  9: savanna
! 10: grassland / steppe
! 11: dense shrubland
! 12: open shrubland
! 13: tundra
! 14: desert
! 15: polar desert / rock / ice
! 16: croplands
! ---------------------------------------------------
!
! these classes consist of some combination of
! plant functional types:
!
! ---------------------------------------------------
!  1: tropical broadleaf evergreen trees
!  2: tropical broadleaf drought-deciduous trees
!  3: warm-temperate broadleaf evergreen trees
!  4: temperate conifer evergreen trees
!  5: temperate broadleaf cold-deciduous trees
!  6: boreal conifer evergreen trees
!  7: boreal broadleaf cold-deciduous trees
!  8: boreal conifer cold-deciduous trees
!  9: evergreen shrubs
! 10: cold-deciduous shrubs
! 11: warm (c4) grasses
! 12: cool (c3) grasses
! 13: soybean
! 14: maize 
! 15: winter and spring wheat 
! ---------------------------------------------------
!*** David Price 2001/05/25. The following code replaces the 450+
!    lines defining plai for each vegtype. Note that values of plai_init are read in as
!    parameters from params.veg. Note also that the declarations
!    of the four local variables plaievgr, plaideci, plaishrub 
!    and plaigrass have been dropped.
               plai(i,1)  = exist(i,1)  / float(ievgr)  * plai_init(1,inveg)
               plai(i,2)  = exist(i,2)  / float(ideci)  * plai_init(2,inveg)
               plai(i,3)  = exist(i,3)  / float(ievgr)  * plai_init(1,inveg)
               plai(i,4)  = exist(i,4)  / float(ievgr)  * plai_init(1,inveg)
               plai(i,5)  = exist(i,5)  / float(ideci)  * plai_init(2,inveg)
               plai(i,6)  = exist(i,6)  / float(ievgr)  * plai_init(1,inveg)
               plai(i,7)  = exist(i,7)  / float(ideci)  * plai_init(2,inveg)
               plai(i,8)  = exist(i,8)  / float(ideci)  * plai_init(2,inveg)
               plai(i,9)  = exist(i,9)  / float(ishrub) * plai_init(3,inveg)
               plai(i,10) = exist(i,10) / float(ishrub) * plai_init(3,inveg)

               if (((inveg.eq.9).or.(inveg.eq.10)).and.(nint(exist(i,11)).eq.1)) then
                     plai(i,11) = exist(i,11) * 0.80 * plai_init(4,inveg)
                     plai(i,12) = exist(i,12) * 0.20 * plai_init(4,inveg)
               else
                  plai(i,11) = exist(i,11) / float(igrass) * plai_init(4,inveg)
                  plai(i,12) = exist(i,12) / float(igrass) * plai_init(4,inveg)
               endif
#ifndef SINGLE_POINT_MODEL
! cropping systems - just set equal to zero
               plai(i,13) = 0.0  ! soybean
               plai(i,14) = 0.0  ! maize
               plai(i,15) = 0.0  ! wheat
               plai(i,16) = 0.0  ! sugarcane
#endif /* SINGLE_POINT_MODEL */
            endif !isimveg.eq.0 .or. isimveg.eq.1

! ************************************************************************
! case (2) assign vegetation characteristics for dynamic vegtation
!          that is initialized with uniform vegetation conditions
! ************************************************************************
!
! specify uniform initial conditions
            if (isimveg.eq.2) then
               plai(i,1)  = exist(i,1)  / float(iupper) * plaiupper
               plai(i,2)  = exist(i,2)  / float(iupper) * plaiupper
               plai(i,3)  = exist(i,3)  / float(iupper) * plaiupper
               plai(i,4)  = exist(i,4)  / float(iupper) * plaiupper
               plai(i,5)  = exist(i,5)  / float(iupper) * plaiupper
               plai(i,6)  = exist(i,6)  / float(iupper) * plaiupper
               plai(i,7)  = exist(i,7)  / float(iupper) * plaiupper
               plai(i,8)  = exist(i,8)  / float(iupper) * plaiupper
               plai(i,9)  = exist(i,9)  / float(ilower) * plailower
               plai(i,10) = exist(i,10) / float(ilower) * plailower
               plai(i,11) = exist(i,11) / float(ilower) * plailower
               plai(i,12) = exist(i,12) / float(ilower) * plailower
#ifndef SINGLE_POINT_MODEL 
               plai(i,13) = exist(i,13) / float(ilower) * plailower
               plai(i,14) = exist(i,14) / float(ilower) * plailower
               plai(i,15) = exist(i,15) / float(ilower) * plailower
               plai(i,16) = exist(i,16) / float(ilower) * plailower
#endif /* SINGLE_POINT_MODEL */
            endif

! ************************************************************************
! for both cases (1) and (2)
! ************************************************************************
!
! set minimum lai for each existing plant type
            plai(i,1)  = max (plai(i,1) , exist(i,1)  * xminlai)
            plai(i,2)  = max (plai(i,2) , exist(i,2)  * xminlai)
            plai(i,3)  = max (plai(i,3) , exist(i,3)  * xminlai)
            plai(i,4)  = max (plai(i,4) , exist(i,4)  * xminlai)
            plai(i,5)  = max (plai(i,5) , exist(i,5)  * xminlai)
            plai(i,6)  = max (plai(i,6) , exist(i,6)  * xminlai)
            plai(i,7)  = max (plai(i,7) , exist(i,7)  * xminlai)
            plai(i,8)  = max (plai(i,8) , exist(i,8)  * xminlai)
            plai(i,9)  = max (plai(i,9) , exist(i,9)  * xminlai)
            plai(i,10) = max (plai(i,10), exist(i,10) * xminlai)
            plai(i,11) = max (plai(i,11), exist(i,11) * xminlai)
            plai(i,12) = max (plai(i,12), exist(i,12) * xminlai)
#ifndef SINGLE_POINT_MODEL 
            plai(i,13) = max (plai(i,13), exist(i,13) * xminlai)
            plai(i,14) = max (plai(i,14), exist(i,14) * xminlai)
            plai(i,15) = max (plai(i,15), exist(i,15) * xminlai)
            plai(i,16) = max (plai(i,16), exist(i,16) * xminlai)
#endif /* SINGLE_POINT_MODEL */

! set sapwood fraction and biomass characteristics
            sapfrac(i) = sapfrac_init ! 0.1 from params.veg
            wood = 0.0
            do 120 j = 1, npft
               cbiol(i,j) = plai(i,j) / specla(i,j)		! Castanho HP, 2013
               cbior(i,j) = 0.5 * cbiol(i,j)
               cbiow(i,j) = 0.0

! crop biomass storage -- stem and grain (fruit)
            if(isimagro .gt. 0)then
               cbios(i,j) = 0.0
               cbiog(i,j) = 0.0
            endif
           
!
               if (j.lt.9) cbiow(i,j) = plai(i,j) * 10.0 / 6.0
#ifdef SINGLE_POINT_MODEL
! FIXME: Explain why this value is forced and hardcoded to the model. 
!       Should this change be used to the 2D model as well? -fzm
               if (j.eq.1) cbiow(i,j) = 18.5 * 5/6
               if (j.eq.2) cbiow(i,j) = 18.5 * 1/6
#endif /* SINGLE_POINT_MODEL */

        if(isimagro .eq. 0) then
               biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
        else
               biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j) + &
                              cbios(i,j) + cbiog(i,j)
        endif
#ifndef SINGLE_POINT_MODEL
             if(isimagro .eq. 0) then
               if (j.le.8) wood = wood + cbiow(i,j)
             else
                  wood = wood + cbiow(i,j)
             endif
#else /* SINGLE_POINT_MODEL */
! FIXME: On 0D version, wood is -always- incremented. Not only for j <= 8.
!       This is potentially an issue to the Agro-Ibis version. - fzm
               wood = wood + cbiow(i,j)
#endif /* SINGLE_POINT_MODEL */
120         continue

! total leaf area for upper and lower canopies
            totlaiu(i)  =  plai(i,1) + plai(i,2) + plai(i,3) + plai(i,4) + &
                           plai(i,5) + plai(i,6) + plai(i,7) + plai(i,8)
#ifndef SINGLE_POINT_MODEL 
            totlail(i)  =  plai(i,9) + plai(i,10) + plai(i,11) + plai(i,12) + &
                           plai(i,13) + plai(i,14) + plai(i,15) + plai(i,16)
#else
            totlail(i)  =  plai(i,9) + plai(i,10) + plai(i,11) + plai(i,12)
#endif /* SINGLE_POINT_MODEL */

            totbiou(i)  = biomass(i,1)+biomass(i,2)+biomass(i,3)+biomass(i,4) + &
                          biomass(i,5)+biomass(i,6)+biomass(i,7)+biomass(i,8)
#ifndef SINGLE_POINT_MODEL 
            totbiol(i)  = biomass(i,9)+biomass(i,10)+biomass(i,11)+biomass(i,12) + &
                          biomass(i,13)+biomass(i,14)+biomass(i,15)+biomass(i,16)
#else
            totbiol(i)  = biomass(i,9)+biomass(i,10)+biomass(i,11)+biomass(i,12)

#endif /* SINGLE_POINT_MODEL */

! fractional cover
#ifndef SINGLE_POINT_MODEL
            fu(i) = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
            fl(i) = totlail(i) / 1.0
#endif /* SINGLE_POINT_MODEL */
            fu(i) = max (dble(0.25), min (dble(0.975), fu(i)))
            fl(i) = max (dble(0.25), min (dble(0.975), fl(i)))

! initial single-sided sai for upper and lower canopies
            sai(i,1)  =  0.050 * totlail(i)
            sai(i,2)  =  0.250 * totlaiu(i)

! initial lai for canopy physics
            lai(i,1) = totlail(i) / fl(i)
            lai(i,2) = totlaiu(i) / fu(i)

! specify canopy height parameters
! calculated as a function of only the vegetative fraction
! of each grid cell

! TODO: it may be desired on 0D version (single point input), to fix these
!     values to a defined (previously observed) value set. For instance, on the
!     2D version, the calculations must be dynamical.
 
          if(croptype(i) .lt. scpft)then

            zbot(i,1) = 0.05
            ztop(i,1) = max (dble(0.25), lai(i,1) * 0.25)
#ifndef SINGLE_POINT_MODEL
            zbot(i,2) = ztop(i,1) + 1.0
            ztop(i,2) = max(zbot(i,2) + 1.00, 2.50 * &
                        totbiou(i) / fu(i) * 0.75)
#else /* SINGLE_POINT_MODEL */
! TODO: As said above, this might be read from somewhere instead of being hard-
!      coded here. - fzm
! Imbuzeiro: Fixed the zbot(i,2) for the LBA_DMIP Protocol
!            zbot(i,2) = 15.0

! Imbuzeiro: Fixed the ztop(i,2) for the LBA_DMIP Protocol
!            ztop(i,2) = 35.0
#endif /* SINGLE_POINT_MODEL */

! constrain ztop to be at least 0.5 meter lower than 
! zbot for upper canopy
            ztop(i,1) = min (ztop(i,1), zbot(i,2) - 0.5)

          else

             zbot(i,1) =  0.01
             ztop(i,1) =  0.50

! constrain ztop to be at least 0.5 meter lower than 
! zbot for upper canopy

            ztop(i,2) =  5.0  !to match the tower, wind
          
         end if !end agro

! ************************************************************************
! case (3) assign vegetation characteristics for model restart
! ************************************************************************
         else    ! else for restart if loop
           if (isimagro .eq. 0)then
            do 130 j = 1, npft
               plai(i,j) = cbiol(i,j) * specla(i,j)				! Castanho HP, 2013
               biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
130         continue

! total leaf area for upper and lower canopies
               totlaiu(i) = plai(i,1) + plai(i,2) + plai(i,3) + plai(i,4) + &
                            plai(i,5) + plai(i,6) + plai(i,7) + plai(i,8)
               totlail(i) = plai(i,9) + plai(i,10) + plai(i,11) + plai(i,12)+ &
                            plai(i,13) + plai(i,14) + plai(i,15) + plai(i,16)
               totbiou(i) = biomass(i,1)+biomass(i,2)+biomass(i,3)+biomass(i,4)+ &
                            biomass(i,5)+biomass(i,6)+biomass(i,7)+biomass(i,8)
               totbiol(i) = biomass(i,9)+biomass(i,10)+biomass(i,11)+biomass(i,12)+ &
                            biomass(i,13)+biomass(i,14)+biomass(i,15)+biomass(i,16)

! initial single-sided sai for upper and lower canopies
               sai(i,1) = 0.050 * totlail(i)
               sai(i,2) = 0.250 * totlaiu(i)

! Lai read from restart file
!
! specify canopy height parameters
! calculated as a function of only the vegetative fraction
! of each grid cell
               zbot(i,1) = 0.05
               ztop(i,1) = max (dble(0.25), lai(i,1) * 0.25)
               zbot(i,2) = 3.0
               ztop(i,2) = max(zbot(i,2) + 1.00, 2.50 * totbiou(i) / fu(i) * 0.75)

! constrain ztop to be at least 0.5 meter lower than 
! zbot for upper canopy
               ztop(i,1) = min (ztop(i,1), zbot(i,2) - 0.5)

             else !cropsums

               zbot(i,1) =  0.01
               ztop(i,1) =  0.50

! constrain ztop to be at least 0.5 meter lower than 
! zbot for upper canopy

               ztop(i,2) =  5.0  !to match the tower, wind
          
           endif !end isimagro

         end if  ! end restart if loop
100   continue

! *********************************************
! assign some physical properties of vegetation
! *********************************************
!
! SANT -for the global crop model, we have to be specified for each culture..
    if(isimagro .gt. 0)then
        chifuz = 0.65  !Pousa - this value was declared to compare with Agro-IBIS
      if (iwheattype .gt. 0) then
        chiflz = 0.65 
      else
        chiflz = -0.2        ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
      endif
! CJK chiflz =  0.0          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
    endif ! check for crop existence

! leaf optical properties were taken from Sellers et al., 1996
! and Bonan, 1995
      oriev(1) = max (-chiflz, dble(0.))
      oriev(2) = max (-chifuz, dble(0.))
      orieh(1) = max ( chiflz, dble(0.))
      orieh(2) = max ( chifuz, dble(0.))

! ***********************
! define rooting profiles
! ***********************
!
! define rooting profiles based upon data published in:
!
! Jackson et al., 1996:  A global analysis of root distributions
! for terrestrial biomes, Oecologia, 108, 389-411.
!
! and
!
! Jackson et al., 1997:  A global budget for fine root biomass,
! surface area, and nutrient contents, Proceedings of the National
! Academy of Sciences, 94, 7362-7366.
!
! rooting profiles are defined by the "beta" parameter
!
! beta1 is assigned to the lower vegetation layer (grasses and shrubs)
! beta2 is assigned to the upper vegetation layer (trees)
!
! according to Jackson et al. (1996, 1997), the values of beta
! typically fall in the following range
!
! note that the 1997 paper specifically discusses the distribution
! of *fine roots* (instead of total root biomass), which may be more
! important for water and nutrient uptake
!
! --------------                 ------------   ------------
! forest systems                 beta2 (1996)   beta2 (1997)
! --------------                 ------------   ------------
! tropical evergreen forest:        0.962          0.972
! tropical deciduous forest:        0.961          0.982
! temperate conifer forest:         0.976          0.980
! temperate broadleaf forest:       0.966          0.967
! all tropical/temperate forest:    0.970
! boreal forest:                    0.943          0.943
! all trees:                                       0.976
!
! -------------------------      ------------   ------------
! grassland / shrub systems      beta1 (1996)   beta1 (1997)
! -------------------------      ------------   ------------
! tropical grassland / savanna:     0.972          0.972
! temperate grassland:              0.943          0.943
! all grasses:                      0.952          0.952
! schlerophyllous shrubs:           0.964          0.950
! all shrubs:                       0.978          0.975
! crops:                            0.961
! desert:                           0.975          0.970
! tundra:                           0.914
!
! --------------                 ------------
! all ecosystems                 beta  (1996)
! --------------                 ------------
! all ecosystems:                   0.966
!
! for global simulations, we typically assign the following
! values to the beta parameters
!
! beta1 = 0.950, which is typical for tropical/temperate grasslands
! beta2 = 0.970, which is typical for tropical/temperate forests
!
! however, these values could be (and should be) further refined
! when using the model for specific regions

! calculate depth in centimeters
      totdepth = 0.0
      do 300 k = 1, nsoilay
         totdepth = totdepth + hsoi(k) * 100.0
300   continue

! Kai modify for total soil depth of any depth totdepth
     dru = min(totdepth, alog(0.003)/dlog(beta2))   ! Kai NOTE dlog double precision natural log
     drl = min(totdepth, alog(0.003)/dlog(beta1))   ! Kai NOTE 
!

! normalization factors
      frootnorm1 = 1. - beta1 ** totdepth
      frootnorm2 = 1. - beta2 ** totdepth

! calculate rooting profiles
      do 400 k = 1, nsoilay
         if (k.eq.1) then
            depth(k) = hsoi(k) * 100.0
            froot(k,1) = 1. - beta1 ** depth(k)
            froot(k,2) = 1. - beta2 ** depth(k)
         else
            depth(k) = depth(k-1) + hsoi(k) * 100.0
            froot(k,1) = (1. - beta1 ** depth(k)) - (1. - beta1 ** depth(k-1))
            froot(k,2) = (1. - beta2 ** depth(k)) - (1. - beta2 ** depth(k-1))
         endif
         froot(k,1) = froot(k,1) / frootnorm1
         froot(k,2) = froot(k,2) / frootnorm2
400   continue
      return
end subroutine iniveg

