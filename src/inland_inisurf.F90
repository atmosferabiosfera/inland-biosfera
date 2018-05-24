#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine inisurf (irestart)
! ---------------------------------------------------------------------
! does initialization for model
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comveg
      use inland_comsoi
      use inland_comatm
      use inland_comhyd
      use inland_combcs
      use inland_comcrop, only:isimagro

      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      integer irestart     ! 0: not restart of INLAND, 1 restart of INLAND

!-----------------------Local variables---------------------------------
      integer i,j            ! loop indice
     
!-----------------------------------------------------------------------
! set physical constants (mks)
      stef  = 5.67e-8 
      vonk  = 0.4
      grav  = 9.80616
      tmelt = 273.16
      hvap  = 2.5104e+6
      hfus  = 0.3336e+6
      hsub  = hvap + hfus
      ch2o  = 4.218e+3
      cice  = 2.106e+3
      cair  = 1.00464e+3
      cvap  = 1.81e+3
      rair  = 287.04
      rvap  = 461.0
      cappa = rair / cair
      rhow  = 1.0e+3

! -----------------------------------------------------------------
! constant atmospheric co2 and o2
! -----------------------------------------------------------------
!      o2conc = 0.209000
!      co2conc = 0.000350

! specify the epsilon value for the model
      epsilon = 1.0e-7
    if(isimagro .gt. 0)then
! wet day / dry day flag initialized to dry day (0)
!
      do 100 i = lbeg, lend 
        iwet(i) = 0
        do 105 j = 1,31
          iwetday(i,j) = 0
          precipday(i,j) = 0
105    continue 
100  continue
    endif ! check for crop existence
! these are already zeroed upon initialization!
      vzero(:)=0.0

! zero flux arrays, and global diagnostic arrays
      asurd(:,:)=0.0
      asuri(:,:)=0.0
      totcondub(:)=0.0
      totconduc(:)=0.0
      totcondls(:)=0.0
      totcondl3(:)=0.0
      totcondl4(:)=0.0
    if (irestart .eq. 0) then
      totcondc3(:)=0.0
      totcondc4(:)=0.0
    endif
      ginvap(:)=0.0
      gsuvap(:)=0.0
      gtrans(:)=0.0
      grunof(:)=0.0
      gdrain(:)=0.0

 

! initialize vegetation prognostic variables
!
! initialize all temperature fields and canopy air to 10 degrees C 
! if not restart
#ifndef SINGLE_POINT_MODEL
#define BASE_TEMPERATURE 283.16
#else /* SINGLE_POINT_MODEL */
#define BASE_TEMPERATURE 298.16
#endif /* SINGLE_POINT_MODEL */
      if (irestart .eq. 0) then
         tu(:)=BASE_TEMPERATURE
         ts(:)=BASE_TEMPERATURE
         tl(:)=BASE_TEMPERATURE
         t12(:)=BASE_TEMPERATURE
         t34(:)=BASE_TEMPERATURE

! initialize weather generator 'memory'
         xstore(:,:)=0.0

! initialize canopy air conditions (used in turvap)
         q12(:)=0.0
         q34(:)=0.0
         tlsub(:)=273.16

! initialize all co2 concentrations (mol/mol)
         ciub(:)=350.0e-06
         ciuc(:)=350.0e-06
         cils(:)=350.0e-06
         cil3(:)=350.0e-06
         cil4(:)=350.0e-06
         cic3(:)=350.0e-06
         cic4(:)=350.0e-06

         csub(:)=350.0e-06
         csuc(:)=350.0e-06
         csls(:)=350.0e-06
         csl3(:)=350.0e-06
         csl4(:)=350.0e-06
         csc3(:)=350.0e-06
         csc4(:)=350.0e-06

! initialize stomatal conductance (mol-h2o/m**2/sec)
         gsub(:)=0.5
         gsuc(:)=0.5
         gsls(:)=0.5
         gsl3(:)=0.5
         gsl4(:)=0.5
         gsc3(:)=0.5
         gsc4(:)=0.5

! initialize soil biogeochemistry variables if not restart
!         clitlm(:)=0.0
!         clitls(:)=0.0
!         clitll(:)=0.0
!         clitrm(:)=0.0
!         clitrs(:)=0.0
!         clitrl(:)=0.0
!         clitwm(:)=0.0
!         clitws(:)=0.0
!         clitwl(:)=0.0
!         totcmic(:)=0.0
!         csoislop(:)=0.0
!         csoislon(:)=0.0
!         csoipas(:)=0.0
      endif
      totlit(:)=0.0
      totnlit(:)=0.0
      totfall(:)=0.0
      totalit(:)=0.0
      totrlit(:)=0.0
      totanlit(:)=0.0
      totrnlit(:)=0.0
      totcsoi(:)=0.0
      totnmic(:)=0.0
      tco2mic(:)=0.0
      tnpptot(:)=0.0
      tneetot(:)=0.0
      tnmin(:)=0.0

! initialize carbon lost to atmosphere due
! to biomass burning
      cdisturb(:)=0.0

! initialize phenology flags if not restart (else, read from restart files)
         greenfracl3(:)=1.0
         greenfracl4(:)=1.0
      if (irestart .eq. 0) then
         tempu(:)=1.0
         templ(:)=1.0
         dropu(:)=1.0
         dropls(:)=1.0
         dropl4(:)=1.0
         dropl3(:)=1.0

! initialize water and snow interception fractions (if not restart (else, read from restart files)
         wliqu(:)=0.0
         wliqs(:)=0.0
         wliql(:)=0.0
         wsnou(:)=0.0
         wsnos(:)=0.0
         wsnol(:)=0.0
      endif

! initialize time-filtered npp to positive value
      su(:)=0.0
      ss(:)=0.0
      sl(:)=0.0

#ifdef SINGLE_POINT_MODEL
! clit** are now read from parameter files.
! FIXME: This have been never documented, why are these parameters initialized
!      like that for 0D model? -fzm
      totcmic(1) = 0.1667
      totcsoi(1) = 14.00
      totalit(1) = 0.4528
      totrlit(1) = 0.0412

      falll(1) = 0.165
      fallr(1) = 0.110
      fallw(1) = 0.227
#endif /* SINGLE_POINT_MODEL */

      return
end subroutine inisurf
