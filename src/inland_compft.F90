#include "inland_config.h"
module inland_compft
!
      implicit none
      public
      save
! ------
! compft
! ------
!
! -------------------------
! plant parameters for inland
! -------------------------

      real*8 tau15,  &     ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
     	     kc15,   &	   ! co2 kinetic parameter (mol/mol)
     	     ko15,   &	   ! o2 kinetic parameter (mol/mol) 
     	     cimax,  &	   ! maximum value for ci (needed for model stability)
     	     alpha3, &	   ! intrinsic quantum efficiency for C3 plants (dimensionless)
     	     theta3, &	   ! photosynthesis coupling coefficient for C3 plants (dimensionless)
             thetac3,&     ! photosynthesis coupling coefficient for C3 crops
     	     beta3,  &	   ! photosynthesis coupling coefficient for C3 plants (dimensionless)
     	     alpha4, &	   ! intrinsic quantum efficiency for C4 plants (dimensionless)
     	     beta4,  &	   ! photosynthesis coupling coefficient for C4 plants (dimensionless)
     	     theta4, &	   ! photosynthesis coupling coefficient for C4 plants (dimensionless)
             thetac4,&     ! photosynthesis coupling coefficient for C4 crops
             betac3, &     ! photosynthesis coupling coefficient for C3 crops 
             betac4        ! photosynthesis coupling coefficient for C4 crops

      real*8, dimension (:,:), allocatable :: vmax_pft   ! Castanho HP, 2013 included npoi dimension
      real*8, dimension (:), allocatable :: vmax_pftp   ! Castanho HP, 2013 keep cte one dimension params in canopy added p in name
      
!       vmax_pft(npoi,npft)  ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]  ! Castanho HP, 2013 included npoi dimension
!       vmax_pftp(npft)  ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) one dimension [not used]  ! Castanho HP, 2013 cte parmas in canopy added p in name

      real*8 gammaub, &    ! leaf respiration coefficient
     	     coefmub, &    ! 'm' coefficient for stomatal conductance relationship
     	     coefbub, &    ! 'b' coefficient for stomatal conductance relationship
     	     gsubmin	   ! absolute minimum stomatal conductance
!
      real*8 gammauc, &    ! leaf respiration coefficient
     	     coefmuc, &    ! 'm' coefficient for stomatal conductance relationship  
     	     coefbuc, &    ! 'b' coefficient for stomatal conductance relationship  
     	     gsucmin	   ! absolute minimum stomatal conductance
!
      real*8 gammals, &    ! leaf respiration coefficient
     	     coefbls, &    ! 'b' coefficient for stomatal conductance relationship 
     	     coefmls, &	   ! 'm' coefficient for stomatal conductance relationship 
     	     gslsmin	   ! absolute minimum stomatal conductance
!
      real*8 gammal3, &    ! leaf respiration coefficient
     	     coefml3, &    ! 'm' coefficient for stomatal conductance relationship 
     	     coefbl3, &    ! 'b' coefficient for stomatal conductance relationship 
     	     gsl3min	   ! absolute minimum stomatal conductance
!
      real*8 gammal4, &    ! leaf respiration coefficient
     	     coefml4, &    ! 'm' coefficient for stomatal conductance relationship
     	     coefbl4, &    ! 'b' coefficient for stomatal conductance relationship
     	     gsl4min	   ! absolute minimum stomatal conductance
!
!      real*8, dimension (:), allocatable :: tauleaf, tauroot, tauwood, tauwood0 ! Castanho HP commented
       
       real*8, dimension (:), allocatable :: tauleaf, tauroot, tauwood0p       !  Castanho HP one dimension parms in canopy added p in name (only for tauwood0)
       real*8, dimension (:,:), allocatable :: tauwood, tauwood0    ! Castanho HP included npoi dimension
     
! Castanho HP, 2013 included dimensions npoi (i) in aleaf, awood, aroot, tauwood, specla, vmax when appropriate bellow

      real*8 gammac3, &    ! leaf respiration coefficient
             coefmc3, &    ! 'm' coefficient for stomatal conductance relationship 
             coefbc3, &    ! 'b' coefficient for stomatal conductance relationship 
             gsc3min       ! absolute minimum stomatal conductance
!
      real*8 gammac4, &    ! leaf respiration coefficient
             coefmc4, &    ! 'm' coefficient for stomatal conductance relationship
             coefbc4, &    ! 'b' coefficient for stomatal conductance relationship
             gsc4min       ! absolute minimum stomatal conductance
!
      real*8, dimension (:), allocatable :: lotemp, hitemp, drought, f1, f2
!      
!       tauleaf(npft),  ! foliar biomass turnover time constant (years)
!       tauroot(npft),  ! fine root biomass turnover time constant (years)
!       tauwood(npoi,npft),  ! wood biomass turnover time constant (years)
!       tauwood0(npoi,npft)  ! normal (unstressed) turnover time for wood biomass (years)
!       tauwood0p(npft)  ! normal (unstressed) turnover time for wood biomass (years) one dimension
!       lotemp(npft),   ! low temperature threshold in tempvm equation
!       hitemp(npft),   ! high temperature threshold in tempvm equation 
!       drought(npft),  ! crop sensitivity to drought parameter 
!       f1(npft),       ! constant used in tempvm equations 
!       f2(npft)        ! constant used in tempvm equations
!
      real*8, dimension (:), allocatable :: TminL, TminU, Twarm, GDD
! 
!       TminL(npft),    ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
!       TminU(npft),    ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
!       Twarm(npft),    ! Temperature of warmest month (lower canopy PFTs)
!       GDD(npft)       ! minimum GDD needed (base 5 C for upper canopy PFTs, 
!
      real*8, dimension (:,:), allocatable :: plai_init
!
!       plai_init(4,15), ! initial total LAI for each vegtype (used in iniveg)
!
      real*8 plaiupper,    &  ! Potental LAI of upper canopy (uniform initial vegetation)
             plailower,    &  ! Potental LAI of lower canopy (uniform initial vegetation)
             xminlai,      &  ! Minimum LAI for each existing PFT
             sapfrac_init, &  ! Initial value of sapwood fraction used for all woody PFTs
             chifuz,       &  ! upper canopy leaf orientation factor
             chiflz,       &  ! lower canopy leaf orientation factor
             beta1,        &  ! parameter for Jackson rooting profile, lower canopy
             beta2            ! parameter for Jackson rooting profile, upper canopy
!
end module inland_compft
