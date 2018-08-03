#include "inland_config.h"
module inland_comsum

      implicit none
      public
      save
!
! ------
! comsum
! ------
!
      integer idoy                ! day of year counter  

      integer, dimension(:), allocatable :: ndtime, nmtime, nytime
!
!     ndtime(numlv)      ! counter for daily average calculations
!     nmtime(numlv)      ! counter for monthly average calculations
!     nytime(numlv)      ! counter for yearly average calculations
!
! Kludge for nytime (real array of length npoi)
!
      real*8, dimension(:), allocatable :: anytime
!
! Area of each gridcell
!

! moved to inland_parameters.F90, as it makes part of compar.h on SAGE.
!     real*8, dimension(:), allocatable :: garea
!
!     garea(npoi)        ! area of each gridcell (m**2)
!
! Total water content (used to check mass balance and equilibrium)
!
      real*8, dimension(:), allocatable :: dwtot, wtotp, wtot
!
!     dwtot(npoi) 
!     wtotp(npoi) 
!     wtot(npoi)     ! total amount of water stored in snow, soil, puddles, and on vegetation (kg_h2o)
!
! daily average fields
!
      real*8, dimension(:), allocatable :: td, adrain, adsnow, adaet,         &
                                           adtrunoff, adsrunoff, addrainage,  &
                                           adsnod, adsnof, adwsoi, adwisoi,   &
                                           adtsoi, adwsoic, adtsoic, adco2mic,&
                                           adco2root, adco2soi, adco2ratio,   &
                                           adnmintot, adtlaysoi, adwlaysoi,   &
                                           adneetot, adtnpptot, adtrans,      &
                                           adevap, adtratio, adwsoi2,         &
                                           adpbio, adpmoi,                    &
                                           adpign, adpfire, adsrate, adabday, &
                                           adburnfrac
                                           !adrh, adud,
!
!     td(npoi)            ! daily average temperature (K)
!     adrain(npoi)        ! daily average rainfall rate (mm/day)
!     adsnow(npoi)        ! daily average snowfall rate (mm/day)
!     adaet(npoi)         ! daily average aet (mm/day)
!     adtrunoff(npoi)     ! daily average total runoff (mm/day)
!     adsrunoff(npoi)     ! daily average surface runoff (mm/day)
!     addrainage(npoi)    ! daily average drainage (mm/day)
!     adsnod(npoi)        ! daily average snow depth (m)
!     adsnof(npoi)        ! daily average snow fraction (fraction)
!     adwsoi(npoi)        ! daily average soil moisture (fraction)
!     adwisoi(npoi)       ! daily average soil ice (fraction)
!     adtsoi(npoi)        ! daily average soil temperature (c)
!     adwsoic(npoi)       ! daily average soil moisture using root profile weighting (fraction)
!     adtsoic(npoi)       ! daily average soil temperature (c) using profile weighting
!     adco2mic(npoi)      ! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
!     adco2root(npoi)     ! daily accumulated co2 respiration from roots (kg_C m-2 /day)
!     adco2soi(npoi)      ! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
!     adco2ratio(npoi)    ! ratio of root to total co2 respiration
!     adnmintot(npoi)     ! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
!     adtlaysoi(npoi)     ! daily average soil temperature (c) of top layer
!     adwlaysoi(npoi)     ! daily average soil moisture of top layer(fraction)
!     adneetot(npoi)      ! daily accumulated net ecosystem exchange of co2 in ecosystem (kg-C/m**2/day)       
!     adtnpptot(npoi)     ! daily accumulated net ecosystem exchange of co2 in ecosystem (kg-C/m**2/day)
!     adpbio               ! daily average fire probability dependence on biomass
!     adpmoi(npoi)         ! daily average fire probability dependence on wetness
!     adpign(npoi)         ! daily average fire probability dependence on ignition
!     adpfire(npoi)        ! daily average fire probability = Pbio*Pmoi*Pign
!     adsrate              ! daily average fire spread rate (km/h)
!     adabday              ! daily average burned area in 1 day (km^2)
!     adburnfrac           ! daily average burned fraction (fraction)
!
!     adtrans(npoi)       ! daily average transpiration (mm/day)
!     adevap(npoi)        ! daily average evaporation (mm/day)
!     adtratio(npoi)      ! daily average transpiration : evaporation ratio 
!     adrh(npoi)          ! daily average rh (percent)
!     adud(npoi)          ! daily average windspeed
!     adwsoi2(npoi)       ! daily average soil moisture (fraction)
      real*8, dimension(:,:), allocatable :: adburnpft
!     adburnpft            ! daily average burned area for each PFT (km^2)

! monthly average fields
!
      real*8, dimension(:), allocatable :: amts2, amrain, amsnow, amqa, amaet, amtransu, amtransl, &
                                           amsuvap, aminvap, amtrunoff, amsrunoff, amdrainage, amwsoi, &
                                           amwisoi, amvwc, amawc, amtsoi, amsnod, amsnof, amlaiu, amlail, &
                                           amsolar, amreflect, amirdown, amirup, amsens, amlatent, amnpptot, &
                                           amneetot, amco2mic, amco2root, amco2soi, amco2ratio, amnmintot,   &
					                       a10tmin, a5tmin, a5td,  amtrans, amtratio, amtotnleach, amno3leach, &
                                           ampbio, ampmoi, ampign, ampfire, amsrate, amabday, amburnfrac
                                           !amalbedo ,amtemp, amcloud, amrh
!
!     amts2(npoi)         ! monthly average 2-m air temperature (K)
!     amrain(npoi)        ! monthly average rainfall rate (mm/day)
!     amsnow(npoi)        ! monthly average snowfall rate (mm/day)
!     amqa(npoi)          ! monthly average specific humidity (kg-h2o/kg-air)
!     amaet(npoi)         ! monthly average aet (mm/day)
!     amtransu(npoi)      ! monthly average transp. of upper canopy (mm/day)
!     amtransl(npoi)      ! monthly average transp. of lower canopy (mm/day)
!     amsuvap(npoi)       ! monthly average evaporation of soil (mm/day)
!     aminvap(npoi)       ! monthly average interception loss (mm/day)
!     amtrunoff(npoi)     ! monthly average total runoff (mm/day)
!     amsrunoff(npoi)     ! monthly average surface runoff (mm/day)
!     amdrainage(npoi)    ! monthly average drainage (mm/day)
!     amwsoi(npoi)        ! monthly average 1m soil moisture (fraction)
!     amwisoi(npoi)       ! monthly average 1m soil ice (fraction)
!     amvwc(npoi)         ! monthly average 1m volumetric water content (fraction)
!     amawc(npoi)         ! monthly average 1m plant-available water content (fraction)
!     amtsoi(npoi)        ! monthly average 1m soil temperature (C)
!     amsnod(npoi)        ! monthly average snow depth (m)
!     amsnof(npoi)        ! monthly average snow fraction (fraction)
!     amlaiu(npoi)        ! monthly average lai for upper canopy (m**2/m**2)
!     amlail(npoi)        ! monthly average lai for lower canopy (m**2/m**2)
!     amsolar(npoi)       ! monthly average incident solar radiation (W/m**2)
!     amreflect(npoi)     ! monthly average reflected solar radiation (W/m**2) 
!     amirdown(npoi)      ! monthly average downward ir radiation (W/m**2)
!     amirup(npoi)        ! monthly average upward ir radiation (W/m**2)
!     amsens(npoi)        ! monthly average sensible heat flux (W/m**2)
!     amlatent(npoi)      ! monthly average latent heat flux (W/m**2)
!     amnpptot(npoi)      ! monthly total npp for ecosystem (kg-C/m**2/month)
!     amneetot(npoi)      ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
!     amco2mic(npoi)      ! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
!     amco2root(npoi)     ! monthly total CO2 flux from soil due to root respiration (kg-C/m**2/month)
!     amco2soi(npoi)      ! monthly total soil CO2 flux from microbial and root respiration (kg-C/m**2/month)
!     amco2ratio(npoi)    ! monthly ratio of root to total co2 flux
!     amnmintot(npoi)     ! monthly total N mineralization from microbes (kg-N/m**2/month)
!     ampbio              ! monthly average fire probability dependence on biomass
!     ampmoi(npoi)        ! monthly average fire probability dependence on wetness
!     ampign(npoi)        ! monthly average fire probability dependence on ignition
!     ampfire(npoi)       ! monthly average fire probability = Pbio*Pmoi*Pign
!     amsrate             ! monthly average fire spread rate (km/h)
!     amabday             ! monthly average burned area in 1 day (km^2)
!     amburnfrac          ! monthly total burned fraction (fraction)
!     a10tmin(npoi)       ! 10-day average minimum air temperature (K)
!     a5tmin(npoi)        ! 5-day average minimum air temperature (K)
!     a5td(npoi)          ! 5-day average daily air temperature (K)
!     amtemp(npoi),       ! monthly average air temperature (C)
!     amcloud(npoi),      ! monthly average cloudiness (percent)
!     amrh(npoi),         ! monthly average rh (percent)
!     amtrans(npoi),      ! monthly average transpiration (mm/day)
!     amtratio(npoi),     ! monthly average transpiration:ET ratio (fraction)
!     amtotnleach(npoi),  ! monthly average total nitrogen leaching (kg/ha/day)
!     amno3leach(npoi),   ! monthly average nitrate leaching (kg/ha/day)
!     amalbedo(npoi),     ! monthly average solar albedo (fraction)
!
      real*8, dimension(:,:), allocatable :: amnpp, aybprod, ayabprod, ayrprod, aylprod, &
                                             adnpp, amburnpft
!
!     amnpp(npoi,npft)    ! monthly total npp for each plant type (kg-C/m**2/month)
!     amburnpft(npoi,npft) ! monthly average burned fraction for each PFT (fraction)
!     aybprod(npoi,npft), ! annual total biomass production for crops
!     ayabprod(npoi,npft),! annual total aboveground biomass production for crop
!     ayrprod(npoi,npft), ! annual total root carbon accumulation for crops (kg-c/m**2/yr) 
!     aylprod(npoi,npft)  ! annual total leaf carbon accumulation for crops (kg-c/m**2/yr)
!     adnpp(npoi,npft)    ! daily total npp for each plant type (kg-C/m**2/day)
!
!
      real*8, dimension(:,:), allocatable :: adcsoln, adwisoilay, adwsoilay
                                             
                                             !adtsoilay, adupsoil
!
!     adcsoln(npoi, nsoilay)     ! daily average nitrate concentration in solution (mg/liter)
!     adwisoilay(npoi,nsoilay),  ! daily average soil ice content for each soil layer
!     adwsoilay(npoi,nsoilay),   ! daily average soil moisture for each soil layer
!     adtsoilay(npoi, nsoilay),  ! daily average soil temperature for each soil layer  
!     adupsoil(npoi, nsoilay),   ! daily total soil water uptake from lower canopy (kg m-2/day)
!
! annual average fields
!
      real*8, dimension(:), allocatable :: ayprcp, ayaet, aytrans, aytrunoff, aytrunoff2, aysrunoff, aydrainage, aywsoi,      &
                                           aywisoi, ayvwc, ayawc, aytsoi, ayrratio, aytratio, aysolar, ayreflect, &
                                           ayirdown, ayirup, aysens, aylatent, aystresstu, aystresstl, ayanpptot, &
                                           aynpptot, aygpptot, ayalit, ayblit, aycsoi, aycmic, ayanlit, aybnlit,  &
                                           aynsoi, ayneetot, ayco2mic, ayco2root, ayco2soi, aynmintot, ayrootbio, &
                                           aymintot, ayimmtot, aynreltot, ayalbedo,                               &
                                           aypbio, aypmoi, aypign, aypfire, aysrate, ayabday, ayburnfrac
!
!     ayprcp(npoi)        ! annual average precipitation (mm/yr)
!     ayaet(npoi)         ! annual average aet (mm/yr)
!     aytrans(npoi)       ! annual average transpiration (mm/yr)
!     aytrunoff(npoi)     ! annual average total runoff (mm/yr)
!     aysrunoff(npoi)     ! annual average surface runoff (mm/yr)
!     aydrainage(npoi)    ! annual average drainage (mm/yr)
!     aywsoi(npoi)        ! annual average 1m soil moisture (fraction)
!     aywisoi(npoi)       ! annual average 1m soil ice (fraction)
!     ayvwc(npoi)         ! annual average 1m volumetric water content (fraction)
!     ayawc(npoi)         ! annual average 1m plant-available water content (fraction)
!     aytsoi(npoi)        ! annual average 1m soil temperature (C)
!     ayrratio(npoi)      ! annual average runoff ratio (fraction)
!     aytratio(npoi)      ! annual average transpiration ratio (fraction)
!     aysolar(npoi)       ! annual average incident solar radiation (w/m**2)
!     ayreflect(npoi)     ! annual average solar albedo (fraction)
!     ayirdown(npoi)      ! annual average downward ir radiation (w/m**2)
!     ayirup(npoi)        ! annual average upward ir radiation (w/m**2)
!     aysens(npoi)        ! annual average sensible heat flux (w/m**2)
!     aylatent(npoi)      ! annual average latent heat flux (w/m**2)
!     aystresstu(npoi)    ! annual average soil moisture stress parameter for upper canopy (dimensionless)
!     aystresstl(npoi)    ! annual average soil moisture stress parameter for lower canopy (dimensionless)
!     ayanpptot(npoi)     ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
!     aynpptot(npoi)      ! annual total npp for ecosystem (kg-c/m**2/yr)
!     aygpptot(npoi)      ! annual total gpp for ecosystem (kg-c/m**2/yr)
!     ayalit(npoi)        ! aboveground litter (kg-c/m**2)
!     ayblit(npoi)        ! belowground litter (kg-c/m**2)
!     aycsoi(npoi)        ! total soil carbon (kg-c/m**2)
!     aycmic(npoi)        ! total soil carbon in microbial biomass (kg-c/m**2)
!     ayanlit(npoi)       ! aboveground litter nitrogen (kg-N/m**2)
!     aybnlit(npoi)       ! belowground litter nitrogen (kg-N/m**2)
!     aynsoi(npoi)        ! total soil nitrogen (kg-N/m**2)
!     ayneetot(npoi)      ! annual total NEE for ecosystem (kg-C/m**2/yr)
!     ayco2mic(npoi)      ! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
!     ayco2root(npoi)     ! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
!     ayco2soi(npoi)      ! annual total soil CO2 flux from microbial and root respiration (kg-C/m**2/yr)
!     aynmintot(npoi)     ! annual total nitrogen mineralization(kg-N/m**2/yr)
!     ayrootbio(npoi)     ! annual average live root biomass (kg-C/m**2/yr)
!     aypbio              ! annual average fire probability dependence on biomass
!     aypmoi(npoi)        ! annual average fire probability dependence on wetness
!     aypign(npoi)        ! annual average fire probability dependence on ignition
!     aypfire(npoi)       ! annual average fire probability = Pbio*Pmoi*Pign
!     aysrate             ! annual average fire spread rate (km/h)
!     ayabday             ! annual average burned area in 1 day (km^2)
!     ayburnfrac          ! annual total burned fraction (fraction)
!     aymintot(npoi),     ! annual total gross nitrogen mineralization (kg-N/m**2/yr) 
!     ayimmtot(npoi),     ! annual total gross nitrogen immobilization (kg-N/m**2/yr) 
!     aynreltot(npoi),    ! annual total non-microbial nitrogen mineral/immobilization (kg-N/m**2/yr) 
!     ayalbedo(npoi),     ! annual average solar albedo (fraction)
!
      real*8, dimension(:,:), allocatable :: ayanpp, aynpp, aygpp, ayburnpft
!
!     ayanpp(npoi,npft)    ! annual above-ground npp for each plant type(kg-c/m**2/yr)
!     aynpp(npoi,npft)     ! annual total npp for each plant type(kg-c/m**2/yr)
!     aygpp(npoi,npft)     ! annual gross npp for each plant type(kg-c/m**2/yr)
!     ayburnpft(npoi,npft) ! annual average burned fraction for each PFT (km^2)
!
! other time average fields
!
      real*8, dimension(:), allocatable :: a10td, a10ancub, a10ancuc, a10ancls, a10ancl3, a10ancl4, &
                                           a10scalparamu, a10scalparaml, a10daylightu, a10daylightl,&
                                           a10ts, a10ancc3, a10ancc4
!
!     a10td(npoi)         ! 10-day average daily air temperature (K)
!     a10ancub(npoi)      ! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
!     a10ancuc(npoi)      ! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
!     a10ancls(npoi)      ! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
!     a10ancl3(npoi)      ! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
!     a10ancl4(npoi)      ! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
!     a10scalparamu(npoi)     ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
!     a10scalparaml(npoi)     ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
!     a10daylightu(npoi)  ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
!     a10daylightl(npoi)  ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
!     a10ts(npoi),        ! 10-day average daily soil temperature in top layer (K)
!     a10ancc3(npoi),     ! 10-day average canopy photosynthesis rate - c3 crops (mol_co2 m-2 s-1)
!     a10ancc4(npoi),     ! 10-day average canopy photosynthesis rate - c4 crops (mol_co2 m-2 s-1)

! biogeochem summations
!
      real*8, dimension(:), allocatable :: storedn, yrleach, ynleach
!
!     storedn(npoi)       ! total storage of N in soil profile (kg_N m-2) 
!     yrleach(npoi)       ! annual total amount C leached from soil profile (kg_C m-2/yr)
!     ynleach(npoi)       ! annual total amount N leached from soil profile (kg_N m-2/yr)
!

end module inland_comsum
