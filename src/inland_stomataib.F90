#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine stomataib(kpti, kptj)
! ----------------------------------- and others to be specified----------------------------------
! TODO: specify only the needed 'use' for every module
      use inland_parameters, only: dtime, epsilon, tempvm_nec, scpft, ecpft
      use inland_control
      use inland_comatm
      use inland_comcrop, only:grnfraccrop, croplive, isimagro
      use inland_comnitr, only:stressn
      use inland_comsatp
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_com1d
      use inland_compft
      use inland_comhour, only:imetyear


      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer i, idc, iter, j
      real*8 rwork, tau, tleaf, tempvm, zweight, scale, extpar
      real*8 esat12, qsat12, rh12, esat34, qsat34, rh34
      real*8 gbco2u, gbco2l, gscub, gscuc, gscls, gscl3, gscl4, gscc3, gscc4
      real*8 vmax, vmaxub, vmaxuc, vmaxls, vmaxl3, vmaxl4
      real*8 rdarkub, rdarkuc, rdarkls, rdarkl3, rdarkl4, rdarkc3, rdarkc4
      real*8 agub, aguc, agls, agl3, agl4, agc3, agc4
      real*8 anub, anuc, anls, anl3, anl4, anc3, anc4
      real*8 duma, dumb, dumc, dume, dumq, dump
      real*8 pxaiu, plaiu, pxail, plail
      real*8 cscub, cscuc, cscls, cscl3, cscl4, cscc3, cscc4
      real*8 stressc3c, stressc4c

      real*8 kc,   &  ! co2 kinetic parameter (mol/mol)
             ko,   &  ! o2  kinetic parameter (mol/mol)
             kco2, &  ! initial c4 co2 efficiency (mol-co2/m**2/s)
             je,   &  ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
             jc,   &  ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
             ji,   &  ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
             jp,   &  ! model-intermediate rate of photosynthesis (mol-co2/m**2/s)
             js,   &  ! sucrose synthesis limitation
             gamstar,&! gamma*, the co2 compensation points for c3 plants
             q10,  &  !
             tkco2

!-----------------------------------------------------------------------
! include water vapor functions
#define STOMATAIB_COMSAT
#include "inland_comsat.h"
! ---------------------------------------------------------------------
! * * * upper canopy physiology calculations * * *
! ---------------------------------------------------------------------
      do 100 i = kpti, kptj

! calculate physiological parameter values which are a function of temperature
         rwork = 3.47e-03 - 1. / tu(i)
         kc  = kc15  * exp( 6000.0 * rwork)
         tleaf = tu(i) - 273.16

         if(isimagro .eq. 0) then
            tau = tau15 * exp(-5000.0 * rwork)
            ko  = ko15  * exp( 1400.0 * rwork)
            tempvm = exp(tempvm_nec * rwork ) / ((1.0+exp(0.40 * (5.0-tleaf))) * &
                        (1.0 + exp(0.40 * (tleaf-50.0))))
         else
            tau = tau15 * exp(-4500.0 * rwork)
            ko  = ko15  * exp( 1500.0 * rwork)
            tempvm = exp(3500.0 * rwork ) / ((1.0+exp(0.40 * (  5.0 - tleaf))) * &
                        (1.0 + exp(0.40 * (tleaf - 50.0))))
         endif
        

! upper canopy gamma-star values (mol/mol)
         gamstar = o2conc / (2. * tau)

! constrain ci values to acceptable bounds -- to help ensure numerical stability
         ciub(i) = max (1.05 * gamstar, min (cimax, ciub(i)))
         ciuc(i) = max (1.05 * gamstar, min (cimax, ciuc(i)))

! calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
         gbco2u = min (dble(10.0), max (dble(0.1), su(i) * 25.5))

! calculate the relative humidity in the canopy air space
! with a minimum value of 0.30 to avoid errors in the 
! physiological calculations
         esat12 = esat (t12(i))
         qsat12 = qsat (esat12, psurf(i))
         rh12   = max (dble(0.30), q12(i) / qsat12)

! ---------------------------------------------------------------------
! broadleaf (evergreen & deciduous) tree physiology 
! ---------------------------------------------------------------------
! 
! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
!
! tropical broadleaf trees          60.0 e-06 mol/m**2/sec
! warm-temperate broadleaf trees    40.0 e-06 mol/m**2/sec
! temperate broadleaf trees         25.0 e-06 mol/m**2/sec
! boreal broadleaf trees            25.0 e-06 mol/m**2/sec

! Castanho HP, 2013 included dimensions npoi (i) in aleaf, awood, aroot, tauwood, specla, vmax when appropriate bellow

         if (exist(i,1).gt.0.5) then
            vmaxub = vmax_pft(i,1) ! 65.0e-06 ! Tropical broadleaf evergreen				! Castanho HP, 2013
         else if (exist(i,3).gt.0.5) then
            vmaxub = vmax_pft(i,3) ! 40.0e-06 ! Warm-temperate broadleaf evergreen			! Castanho HP, 2013
         else 
            vmaxub = vmax_pft(i,5) ! 30.0e-06 ! Temperate or boreal broadleaf cold deciduous	! Castanho HP, 2013
         endif

! vmax and dark respiration for current conditions
         vmax  = vmaxub * tempvm * stresstu(i)
         rdarkub = gammaub * vmaxub * tempvm

! 'light limited' rate of photosynthesis (mol/m**2/s)
         je = topparu(i) * 4.59e-06 * alpha3 * (ciub(i) - gamstar) / (ciub(i) + 2. * gamstar)

! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
         jc = vmax * (ciub(i) - gamstar) / (ciub(i) + kc * (1. + o2conc / ko))

! solution to quadratic equation
         duma = theta3
         dumb = je + jc
         dumc = je * jc
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the intermediate photosynthesis rate (mol/m**2/s)
         jp = min (dumq/duma, dumc/dumq)

! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
         js = vmax / 2.2

! solution to quadratic equation
         duma = beta3
         dumb = jp + js
         dumc = jp * js
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the net photosynthesis rate (mol/m**2/s)
         agub = min (dumq/duma, dumc/dumq)
         anub = agub - rdarkub

! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
         csub(i) = 0.5 * (csub(i) + co2conc - anub / gbco2u)
         csub(i) = max (1.05 * gamstar, csub(i))

! calculate new value of gs using implicit scheme
         gsub(i) = 0.5 * (gsub(i) + (coefmub * anub * rh12 / csub(i) + &
                                     coefbub * stresstu(i)))
         gsub(i) = max (gsubmin, coefbub * stresstu(i), gsub(i))

! calculate new value of ci using implicit scheme
         ciub(i) = 0.5 * (ciub(i) + csub(i) - 1.6 * anub / gsub(i))
         ciub(i) = max (1.05 * gamstar, min (cimax, ciub(i)))

! ---------------------------------------------------------------------
! conifer tree physiology 
! ---------------------------------------------------------------------
!
! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
!
! temperate conifer trees           30.0 e-06 mol/m**2/sec
! boreal conifer trees              20.0 e-06 mol/m**2/sec
         if (exist(i,4).gt.0.5) then
            vmaxuc = vmax_pft(i,4) ! 30.0e-06 ! Temperate conifer		! Castanho HP, 2013
         else 
            vmaxuc = vmax_pft(i,6) ! 25.0e-06 ! Boreal conifer evergreen	! Castanho HP, 2013
         endif

! vmax and dark respiration for current conditions
         vmax  = vmaxuc * tempvm * stresstu(i)
         rdarkuc = gammauc * vmaxuc * tempvm

! 'light limited' rate of photosynthesis (mol/m**2/s)
         je = topparu(i) * 4.59e-06 * alpha3 * (ciuc(i) - gamstar) / &
             (ciuc(i) + 2. * gamstar)

! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
         jc = vmax * (ciuc(i) - gamstar) / (ciuc(i) + kc * (1. + o2conc / ko))

! solution to quadratic equation
         duma = theta3
         dumb = je + jc
         dumc = je * jc
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the intermediate photosynthesis rate (mol/m**2/s)
         jp = min (dumq/duma, dumc/dumq)

! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
         js = vmax / 2.2

! solution to quadratic equation
         duma = beta3
         dumb = jp + js
         dumc = jp * js
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the net photosynthesis rate (mol/m**2/s)
         aguc = min (dumq/duma, dumc/dumq) 
         anuc = aguc - rdarkuc

! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
         csuc(i) = 0.5 * (csuc(i) + co2conc - anuc / gbco2u)
         csuc(i) = max (1.05 * gamstar, csuc(i))

! calculate new value of gs using implicit scheme
         gsuc(i) = 0.5 * (gsuc(i) + (coefmuc * anuc * rh12 / csuc(i) + &
                                     coefbuc * stresstu(i)))
         gsuc(i) = max (gsucmin, coefbuc * stresstu(i), gsuc(i))

! calculate new value of ci using implicit scheme
         ciuc(i) = 0.5 * (ciuc(i) + csuc(i) - 1.6 * anuc / gsuc(i))
         ciuc(i) = max (1.05 * gamstar, min (cimax, ciuc(i)))

! ---------------------------------------------------------------------
! upper canopy scaling
! ---------------------------------------------------------------------
!
! the canopy scaling algorithm assumes that the net photosynthesis
! is proportional to absored par (apar) during the daytime. during night,
! the respiration is scaled using a 10-day running-average daytime canopy
! scaling parameter.
!
! apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
! an(x) is proportional to apar(x)
!
! therefore, an(x) = an(0) * apar(x) / apar(0)
! an(x) = an(0) * (A exp(-k x) + B exp(-h x) + C exp(h x)) / 
!                 (A + B + C)
!
! this equation is further simplified to
! an(x) = an(0) * exp (-extpar * x)
!
! an(0) is calculated for a sunlit leaf at the top of the canopy using
! the full-blown plant physiology model (Farquhar/Ball&Berry, Collatz).
! then the approximate par extinction coefficient (extpar) is calculated
! using parameters obtained from the two-stream radiation calculation.
!
! an,canopy avg.= integral (an(x), from 0 to xai) / lai
!               = an(0) * (1 - exp (-extpar * xai )) / (extpar * lai)
!
! the term '(1 - exp (-extpar * xai )) / lai)' scales photosynthesis from leaf
! to canopy level (canopy average) at day time. A 10-day running mean of this
! scaling parameter (weighted by light) is then used to scale the respiration
! during night time.
!
! once canopy average photosynthesis is calculated, then the canopy average
! stomatal conductance is calculated using the 'big leaf approach',i.e. 
! assuming that the canopy is a big leaf and applying the leaf-level stomatal
! conductance equations to the whole canopy.
!
! calculate the approximate par extinction coefficient:
!
! extpar = (k * A + h * B - h * C) / (A + B + C)
         extpar = (termu(i,6) * scalcoefu(i,1) + termu(i,7) * scalcoefu(i,2) - &
                   termu(i,7) * scalcoefu(i,3)) / max (scalcoefu(i,4), epsilon)
         extpar = max (dble(1.e-1), min (dble(1.e+1), extpar))

! calculate canopy average photosynthesis (per unit leaf area):
         pxaiu = extpar * (lai(i,2) + sai(i,2))
         plaiu = extpar *  lai(i,2)

! scale is the parameter that scales from leaf-level photosynthesis to
! canopy average photosynthesis
         zweight = exp(-1. / (10.0 * 86400./dtime))

! for non-zero lai
         if (plaiu.gt.0.0) then

! day-time conditions, use current scaling coefficient
            if (topparu(i).gt.10.) then
               scale = (1. - exp(-pxaiu)) / plaiu

! update 10-day running mean of scale, weighted by light levels
               a10scalparamu(i) = zweight * a10scalparamu(i) + &
                                  (1. - zweight) * scale * topparu(i)
               a10daylightu(i) = zweight * a10daylightu(i) + &
                                 (1. - zweight) * topparu(i)

! night-time conditions, use long-term day-time average scaling coefficient
            else
               scale = a10scalparamu(i) / a10daylightu(i)
            endif

! if no lai present
         else
            scale = 0.0
         endif

! perform scaling on all carbon fluxes from upper canopy
         agcub(i) = agub * scale
         agcuc(i) = aguc * scale
         ancub(i) = anub * scale
         ancuc(i) = anuc * scale

! calculate diagnostic canopy average surface co2 concentration
! (big leaf approach)
         cscub = max (1.05 * gamstar, co2conc - ancub(i) / gbco2u)
         cscuc = max (1.05 * gamstar, co2conc - ancuc(i) / gbco2u)

! calculate diagnostic canopy average stomatal conductance (big leaf approach)
         gscub = coefmub * ancub(i) * rh12 / cscub + coefbub * stresstu(i)
         gscuc = coefmuc * ancuc(i) * rh12 / cscuc + coefbuc * stresstu(i)
         gscub = max (gsubmin, coefbub * stresstu(i), gscub)
         gscuc = max (gsucmin, coefbuc * stresstu(i), gscuc)

! calculate total canopy and boundary-layer total conductance for 
! water vapor diffusion
         rwork = 1. / su(i)
         dump  = 1. / 0.029
         totcondub(i) = 1. / (rwork + dump / gscub)
         totconduc(i) = 1. / (rwork + dump / gscuc)

! multiply canopy photosynthesis by wet fraction - this calculation is
! done here and not earlier to avoid using within canopy conductance
         rwork = 1 - fwetu(i)
         agcub(i) = rwork * agcub(i)
         agcuc(i) = rwork * agcuc(i)
         ancub(i) = rwork * ancub(i)
         ancuc(i) = rwork * ancuc(i)
100   continue

! ---------------------------------------------------------------------
! * * * lower canopy physiology calculations * * * *
! ---------------------------------------------------------------------
      do 200 i = kpti, kptj

! calculate physiological parameter values which are a function of temperature
         rwork = 3.47e-03 - 1. / tl(i)
         if(isimagro .eq. 0) then
             tau = tau15 * exp(-4500.0 * rwork)
             ko  = ko15  * exp( 1500.0 * rwork)
         else
             tau = tau15 * exp(-5000.0 * rwork)
             ko  = ko15  * exp( 1400.0 * rwork)
         endif

         kc  = kc15  * exp( 6000.0 * rwork)
         tleaf = tl(i) - 273.16
         tempvm = exp(3500.0 * rwork ) / &
                     ((1.0 + exp(0.40 * (  5.0 - tleaf))) * &
                     (1.0 + exp(0.40 * (tleaf - 50.0))))

! lower canopy gamma-star values (mol/mol)
         gamstar = o2conc / (2. * tau)

! constrain ci values to acceptable bounds -- to help ensure numerical stability
         cils(i) = max (1.05 * gamstar, min (cimax, cils(i)))
         cil3(i) = max (1.05 * gamstar, min (cimax, cil3(i)))
         cil4(i) = max (dble(0.0)           , min (cimax, cil4(i)))
!
! constrain ci value to acceptable bounds for crops
!
       if(isimagro .gt. 0) then
        cic3(i) = max (1.05 * gamstar, min (cimax, cic3(i)))
        cic4(i) = max (0.0           , min (cimax, cic4(i)))
       endif

! calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
         gbco2l = min (dble(10.0), max (dble(0.1), sl(i) * 25.5))

! calculate the relative humidity in the canopy air space
! with a minimum value of 0.30 to avoid errors in the 
! physiological calculations
         esat34 = esat (t34(i))
         qsat34 = qsat (esat34, psurf(i))
         rh34   = max (dble(0.30), q34(i) / qsat34)

! ---------------------------------------------------------------------
! shrub physiology
! ---------------------------------------------------------------------
! 
! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
         vmaxls = vmax_pft(i,9) ! 27.5e-06 ! Shrubs (evergreen or cold deciduous) 		! Castanho HP, 2013

! vmax and dark respiration for current conditions
         vmax  = vmaxls * tempvm * stresstl(i)
         rdarkls = gammals * vmaxls * tempvm

! 'light limited' rate of photosynthesis (mol/m**2/s)
         je = topparl(i) * 4.59e-06 * alpha3 * (cils(i) - gamstar) / &
              (cils(i) + 2. * gamstar)

! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
         jc = vmax * (cils(i) - gamstar) / (cils(i) + kc * (1. + o2conc / ko))

! solution to quadratic equation
         duma = theta3
         dumb = je + jc
         dumc = je * jc
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the intermediate photosynthesis rate (mol/m**2/s)
         jp = min (dumq/duma, dumc/dumq)

! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
         js = vmax / 2.2

! solution to quadratic equation
         duma = beta3
         dumb = jp + js
         dumc = jp * js
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the net photosynthesis rate (mol/m**2/s)
         agls = min (dumq/duma, dumc/dumq)
         anls = agls - rdarkls

! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
         csls(i) = 0.5 * (csls(i) + co2conc - anls / gbco2l)
         csls(i) = max (1.05 * gamstar, csls(i))

! calculate new value of gs using implicit scheme
         gsls(i) = 0.5 * (gsls(i) + coefmls * anls * rh34 / csls(i) + &
                                    coefbls * stresstl(i))
         gsls(i) = max (gslsmin, coefbls * stresstl(i), gsls(i))

! calculate new value of ci using implicit scheme
         cils(i) = 0.5 * (cils(i) + csls(i) - 1.6 * anls / gsls(i))
         cils(i) = max (1.05 * gamstar, min (cimax, cils(i)))

! ---------------------------------------------------------------------
! c3 grass physiology
! ---------------------------------------------------------------------
! 
! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
         vmaxl3 = vmax_pft(i,12) ! 25.0e-06 ! C3 grasses	! 

! vmax and dark respiration for current conditions
         vmax  = vmaxl3 * tempvm * stresstl(i)
         rdarkl3 = gammal3 * vmaxl3 * tempvm

! 'light limited' rate of photosynthesis (mol/m**2/s)
         je = topparl(i) * 4.59e-06 * alpha3 * (cil3(i) - gamstar) / &
              (cil3(i) + 2. * gamstar)

! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
         jc = vmax * (cil3(i) - gamstar) / (cil3(i) + kc * (1. + o2conc / ko))

! solution to quadratic equation
         duma = theta3
         dumb = je + jc
         dumc = je * jc
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the intermediate photosynthesis rate (mol/m**2/s)
         jp = min (dumq/duma, dumc/dumq)

! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
         js = vmax / 2.2

! solution to quadratic equation
         duma = beta3
         dumb = jp + js
         dumc = jp * js
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the net photosynthesis rate (mol/m**2/s)
         agl3 = min (dumq/duma, dumc/dumq)
         anl3 = agl3 - rdarkl3

! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
         csl3(i) = 0.5 * (csl3(i) + co2conc - anl3 / gbco2l)
         csl3(i) = max (1.05 * gamstar, csl3(i))

! calculate new value of gs using implicit scheme
         gsl3(i) = 0.5 * (gsl3(i) + coefml3 * anl3 * rh34 / csl3(i) + &
                          coefbl3 * stresstl(i))
         gsl3(i) = max (gsl3min, coefbl3 * stresstl(i), gsl3(i))

! calculate new value of ci using implicit scheme
         cil3(i) = 0.5 * (cil3(i) + csl3(i) - 1.6 * anl3 / gsl3(i))
         cil3(i) = max (1.05 * gamstar, min (cimax, cil3(i)))

! ---------------------------------------------------------------------
! c4 grass physiology
! ---------------------------------------------------------------------
!
! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
         vmaxl4 = vmax_pft(i,11) ! 15.0e-06 ! C4 grasses		! Castanho HP, 2013

! calculate the parameter values which are a function of temperature
         rwork = 3.47e-03 - 1. / tl(i)
         tleaf = tl(i) - 273.16
         tempvm = exp(3500.0 * rwork ) / &
                  ((1.0 + exp(0.40 * ( 10.0 - tleaf))) * &
                   (1.0 + exp(0.40 * (tleaf - 50.0))))

! vmax and dark respiration for current conditions
         vmax  = vmaxl4 * tempvm * stresstl(i)
         rdarkl4 = gammal4 * vmaxl4 * tempvm

! initial c4 co2 efficiency (mol/m**2/s)
         kco2 = 18.0e+03 * vmax

! 'light limited' rate of photosynthesis (mol/m**2/s)
         je = topparl(i) * 4.59e-06 * alpha4

! 'rubisco limited' rate of photosynthesis
         jc = vmax

! solve for intermediate photosynthesis rate
         duma = theta4
         dumb = je + jc
         dumc = je * jc
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
         jp = min (dumq/duma, dumc/dumq)

! 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
         ji = kco2 * cil4(i)

! solution to quadratic equation
         duma = beta4
         dumb = jp + ji
         dumc = jp * ji
         dume = max (dumb**2 - dble(4.) * duma * dumc, dble(0.))
         dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

! calculate the net photosynthesis rate (mol/m**2/s)
         agl4 = min (dumq/duma, dumc/dumq)
         anl4 = agl4 - rdarkl4

! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
         csl4(i) = 0.5 * (csl4(i) + co2conc - anl4 / gbco2l)
         csl4(i) = max (dble(1.0e-8), csl4(i))

! calculate new value of gs using implicit scheme
         gsl4(i) = 0.5 * (gsl4(i) + coefml4 * anl4 * rh34 / csl4(i) + &
                          coefbl4 * stresstl(i))
         gsl4(i) = max (gsl4min, coefbl4 * stresstl(i), gsl4(i))

! calculate new value of ci using implicit scheme
         cil4(i) = 0.5 * (cil4(i) + csl4(i) - 1.6 * anl4 / gsl4(i))
         cil4(i) = max (dble(0.0), min (cimax, cil4(i)))
    if(isimagro .gt. 0) then
!
!
! get index for current cropping practice - could have more than
! one crop existing in grid cell during the year for multiple
! cropping - but only one would be in live vegetation stage 
!
       idc = 0
       do 250 j = scpft, ecpft 
         if (exist(i,j) .eq. 1. .and. croplive(i,j) .gt. 0) then
            idc = j 
         endif
 250   continue 
!
! --------------------------------------------------------------------
! c3 crops physiology (soybean, wheat)
! ---------------------------------------------------------------------
       if (idc .eq. 13 .or. idc .eq. 15) then  ! soybean or wheat  
!
         rwork = 3.47e-03 - 1. / tl(i)
!
         tleaf = tl(i) - 273.16
!
         q10 = 2.6
!
! vmax and dark respiration for current conditions
!
!         tempvm = exp(3500.0 * rwork ) /
!     >            ((1.0 + exp(0.40 * (  lotemp(idc) - tleaf))) * 
!     >            (1.0 + exp(0.40 * (tleaf - hitemp(idc)))))
!
        tempvm = q10**((tleaf-15.0)/10.0) /   & ! Collatz approach
                 ((1.0 + exp(f1(idc) * (lotemp(idc) - tleaf))) * &
                 (1.0 + exp(f2(idc) * (tleaf - hitemp(idc)))))

!
! adjust drystress factor imposed on soybeans - on a scale of
! 0 (less) - 1.0 (more), these have a 0.8 rating compared to 0.65 for maize 
! and for wheat
! from Penning de Vries, "Simulation of ecophysiological processes of
! growth in several annual crops"
!
!       NEW CJK 9-24-04
        stressc3c = 1.0
	stressn(i,idc)= 1.0
        vmax      = max(0., vmax_pftp(idc) * tempvm * &
                    min(stressc3c, stressn(i,idc), croplive(i,idc)))

        rdarkc3   = gammac3 * vmax_pftp(idc) * tempvm * croplive(i,idc)
!
! 'light limited' rate of photosynthesis (mol/m**2/s)
!
        je = topparl(i) * 4.59e-06 * alpha3 * (cic3(i) - gamstar) / &
             (cic3(i) + 2. * gamstar)
!
! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
!
        jc = vmax * (cic3(i) - gamstar) / &
             (cic3(i) + kc * (1. + o2conc / ko))
!
        duma = thetac3
        dumb = je + jc
        dumc = je * jc
!
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
!
! calculate the intermediate photosynthesis rate (mol/m**2/s)
!
        jp = min (dumq/duma, dumc/dumq)
!       
! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
!       
        js = vmax / 2.2
! 
! solution to quadratic equation
!
        duma = betac3
        dumb = jp + js
        dumc = jp * js
!       
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
!
! calculate the net photosynthesis rate (mol/m**2/s)
!
        agc3 = min (dumq/duma, dumc/dumq)
        anc3 = agc3 - rdarkc3

!
! apply stress functions to net photosynthesis rate
!
        anc3 = anc3 * stresstl(i)  ! CJK 6/20/2004 

!
! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
!
        csc3(i) = 0.5 * (csc3(i) + co2conc - anc3 / gbco2l)
        csc3(i) = max (1.05 * gamstar, csc3(i))
!
! calculate new value of gs using implicit scheme
!
        gsc3(i) = 0.5 * (gsc3(i) + coefmc3 * anc3 * rh34 / csc3(i) + &
                         coefbc3 * stressc3c)
!
        gsc3(i) = max (gsc3min, coefbc3 * stressc3c, gsc3(i))
!
! calculate new value of ci using implicit scheme
!
        cic3(i) = 0.5 * (cic3(i) + csc3(i) - 1.6 * anc3 / gsc3(i))
        cic3(i) = max (1.05 * gamstar, min (cimax, cic3(i)))
!
      else
        agc3    = 0.
        anc3    = 0.     
        csc3(i) = 0.
        gsc3(i) = 0.
        cic3(i) = 0.
!
      endif  ! c3 crops

!      if(dtime.eq.3600*12.and.i.eq.1.and.iter.eq.3.and.iyear.eq.2010.and.jday.le.120) then
      if(imetyear .ne. 9999 .and. idc .ne. 0) then
              write(23,134),jday,dtime/3600.,lai(i,1)*fl(i)*grnfraccrop(i,idc),               &
                            js*1e+06,je*1e+06, jc*1e+06, anc3*1e+06,tleaf,tempvm
      endif
134   format (1x,i3,1x,11(1x,f6.2))	
   
!
! ---------------------------------------------------------------------
! c4 crop physiology (corn and sugarnace)
!
! modified by CJK to follow Collatz et al. (1992) parameterizations
! that were based on corn data in the field
! ---------------------------------------------------------------------
!
      if (idc .eq. 14 .or. idc .eq. 16 ) then   ! maize or sugarcane    
!
! calculate the parameter values which are a function of temperature
!
!sant-original Maize         q10   = 2.0
              q10   = 2.5
        rwork = 3.47e-03 - 1. / tl(i)
        tleaf = tl(i) - 273.16
!
        tempvm = q10**((tleaf-15.0)/10.0) /  &  ! Collatz approach
                 ((1.0 + exp(f1(idc) * (lotemp(idc) - tleaf))) * &
                 (1.0 + exp(f2(idc) * (tleaf - hitemp(idc)))))
!
        stressc4c = 1.0   ! CJK 6/20/2004

 	if(idc.eq.16) stressn(i,idc)=1.0  

        vmax    = max(0., vmax_pftp(idc) * tempvm * &
                  min(stressc4c, stressn(i,idc), croplive(i,idc)))
!
        rdarkc4 = gammac4 * vmax_pftp(idc) * tempvm * croplive(i,idc)

        kco2 = 4.0e+03  * vmax

!
! 'light limited' rate of photosynthesis (mol/m**2/s)
!
        je = topparl(i) * 4.59e-06  * 0.067 ! needed to increase efficiency of corn plants 
!
! 'rubisco limited' rate of photosynthesis
!
        jc = vmax
!
! solve for intermediate photosynthesis rate
!
        duma = thetac4
        dumb = je + jc
        dumc = je * jc
!
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
!
        jp = min (dumq/duma, dumc/dumq)
!
! 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
!
        ji = kco2 * cic4(i)
!
        duma = betac4
        dumb = jp + ji
        dumc = jp * ji
!
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
!
! calculate the net photosynthesis rate (mol/m**2/s)
!
        agc4 = min (dumq/duma, dumc/dumq)
        anc4 = agc4 - rdarkc4
!
! apply stress functions to net photosynthesis rate

        anc4 = anc4 * max(0.0, stresstl(i))  ! CJK 6/20/2004
!
! calculate co2 concentrations and stomatal condutance values
! using simple iterative procedure
!
! weight results with the previous iteration's values -- this
! improves convergence by avoiding flip-flop between diffusion
! into and out of the stomatal cavities
!
! calculate new value of cs using implicit scheme
!

         csc4(i) = (1./3.) * ( csc4(i)*(3.-1.) + co2conc - anc4 / gbco2l)
         csc4(i) = max (0.0, csc4(i))

!
! calculate new value of gs using implicit scheme
!


        gsc4(i) = 0.5 * (gsc4(i) + coefmc4 * anc4 * rh34 / csc4(i) + &
                        coefbc4 * stressc4c)
!

        gsc4(i) = max (gsc4min, coefbc4 * stressc4c, gsc4(i)) 

!
! calculate new value of ci using implicit scheme
!

        cic4(i) = (1/3.) * (cic4(i)*2. + csc4(i) - 1.6 * anc4 / gsc4(i))

        cic4(i) = max (0.0, min (cimax, cic4(i)))
!
      else
        agc4    = 0.
        anc4    = 0.     
        csc4(i) = 0.
        gsc4(i) = 0.
        cic4(i) = 0.
!
      endif    ! c4 crops

    endif ! check for crop existence

! ---------------------------------------------------------------------
! lower canopy scaling
! ---------------------------------------------------------------------
!
! calculate the approximate extinction coefficient
         extpar = (terml(i,6) * scalcoefl(i,1) + &
                  terml(i,7) * scalcoefl(i,2) -  &
                  terml(i,7) * scalcoefl(i,3)) / &
                  max (scalcoefl(i,4), epsilon)
         extpar = max (dble(1.e-1), min (dble(1.e+1), extpar))

! calculate canopy average photosynthesis (per unit leaf area):
         pxail = extpar * (lai(i,1) + sai(i,1))
         plail = extpar *  lai(i,1)

! scale is the parameter that scales from leaf-level photosynthesis to
! canopy average photosynthesis
         zweight = exp(-1. / (10.0 * 86400./dtime))

! for non-zero lai
         if (plail.gt.0.0) then

! day-time conditions, use current scaling coefficient
            if (topparl(i).gt.10.) then
               scale = (1. - exp(-pxail)) / plail

! update 10-day running mean of scale, weighted by light levels
               a10scalparaml(i) = zweight * a10scalparaml(i) + &
                                  (1. - zweight) * scale * topparl(i)
               a10daylightl(i)  = zweight * a10daylightl(i) + &
                                  (1. - zweight) * topparl(i)

! night-time conditions, use long-term day-time average scaling coefficient
            else
               scale = a10scalparaml(i) / a10daylightl(i)
            endif

! if no lai present
         else
            scale = 0.0
         endif

! perform scaling on all carbon fluxes from upper canopy
         agcls(i) = agls * scale
         agcl3(i) = agl3 * scale
       if(isimagro .gt. 0) then
         agcc3(i) = agc3 * scale
         agcc4(i) = agc4 * scale
       endif
         agcl4(i) = agl4 * scale
         ancls(i) = anls * scale
         ancl3(i) = anl3 * scale
         ancl4(i) = anl4 * scale
       if(isimagro .gt. 0) then
         ancc3(i) = anc3 * scale
         ancc4(i) = anc4 * scale
       endif

! calculate canopy average surface co2 concentration
         cscls = max (1.05 * gamstar, co2conc - ancls(i) / gbco2l)
         cscl3 = max (1.05 * gamstar, co2conc - ancl3(i) / gbco2l)
         cscl4 = max (dble(1.0e-8)        , co2conc - ancl4(i) / gbco2l)
!
       if(isimagro .gt. 0) then
         if (idc .eq. 13 .or. idc .eq. 15) then
            cscc3 = max (1.05 * gamstar, co2conc - ancc3(i) / gbco2l)
            cscc4 = 0.
         else if (idc .eq. 14 .or. idc.eq.16) then
            cscc4 = max (1.0e-08       , co2conc - ancc4(i) / gbco2l)
            cscc3 = 0.
         endif
       endif 

! calculate canopy average stomatal conductance
         gscls = coefmls * ancls(i) * rh34 / cscls + coefbls * stresstl(i)
         gscl3 = coefml3 * ancl3(i) * rh34 / cscl3 + coefbl3 * stresstl(i)
         gscl4 = coefml4 * ancl4(i) * rh34 / cscl4 + coefbl4 * stresstl(i)
      if(isimagro .gt. 0) then
!
! c3/c4 crop physiology
!
        if (idc .eq. 13 .or. idc .eq. 15) then
          gscc3 = coefmc3 * ancc3(i) * rh34 / cscc3 + &
                  coefbc3 * stressc3c
          gscc4 = 0.
!
        else if (idc .eq. 14 .or. idc.eq.16) then
          gscc4 = coefmc4 * ancc4(i) * rh34 / cscc4 + &
                  coefbc4 * stressc4c
          gscc3 = 0.  
        else
          gscc3 = 0.
          gscc4 = 0.
        endif

      endif ! check for crop existence

         gscls = max (gslsmin, coefbls * stresstl(i), gscls)
         gscl3 = max (gsl3min, coefbl3 * stresstl(i), gscl3)
         gscl4 = max (gsl4min, coefbl4 * stresstl(i), gscl4)

      if(isimagro .gt. 0) then
!
! c3/c4 crop physiology
!
        if (idc .eq. 13 .or. idc .eq. 15) then
          gscc3 = max (gsc3min, coefbc3 * stressc3c, gscc3)
          gscc3 = gscc3 * grnfraccrop(i,idc)
          gscc4 = 0.
        else if (idc.eq.14 .or. idc.eq.16) then
        ancc4(i) = ancc4(i) * grnfraccrop(i,idc)
        agcc4(i) = agcc4(i) * grnfraccrop(i,idc)

          gscc4 = max (gsc4min, coefbc4 * stressc4c, gscc4)
          gscc4 = gscc4 * grnfraccrop(i,idc)
          gscc3 = 0.
        else
          gscc3 = 0.
          gscc4 = 0.
        endif

!
! The following adjusts the above calculated values of ancl3, ancl4,
! agcl3, agcl4, gscl3, and gscl4 according to what percentage of the
! lower canopy is green by weighting the above calculations by greenfrac
! terms. Only the green portion of the canopy performs photosynthesis.
! Shrubs that have leaves have only green leaves since they are allowed
! to drop their leaves in the fall. C3 and C4 grasses may be either green
! or brown so they are the only terms that are adjusted.
!
! Scale value of ancl3, ancl4, gscl3, and gscl4 according to what fraction
! of the canopy is green
!
        ancl3(i) = ancl3(i) * greenfracl3(i)
        ancl4(i) = ancl4(i) * greenfracl4(i)
!
        agcl3(i) = agcl3(i) * greenfracl3(i)
        agcl4(i) = agcl4(i) * greenfracl4(i)
!
        gscl3 = gscl3 * greenfracl3(i)
        gscl4 = gscl4 * greenfracl4(i)
    endif ! check for crop existence
!
! calculate canopy and boundary-layer total conductance for water vapor diffusion
         rwork = 1. / sl(i)
         dump =  1. / 0.029
!
        totcondls(i) = 1. / (rwork + dump / gscls)
!
! Make sure that the calculation does not divide by zero if gscl3 or
! gscl4 are equal to zero
!
        if (gscl3 .gt. 0) then
          totcondl3(i) = 1. / ( rwork + dump / gscl3 )
        else
          totcondl3(i) = 0
        endif
!
        if (gscl4 .gt. 0) then
          totcondl4(i) = 1. / ( rwork + dump / gscl4 )
        else
          totcondl4(i) = 0
        endif
!
        if (gscc3 .gt. 0) then
          totcondc3(i) = 1. / ( rwork + dump / gscc3 )
        else
          totcondc3(i) = 0
        endif
!
        if (gscc4 .gt. 0) then
          totcondc4(i) = 1. / ( rwork + dump / gscc4 )
        else
          totcondc4(i) = 0
        endif

! multiply canopy photosynthesis by wet fraction -- this calculation is
! done here and not earlier to avoid using within canopy conductance
         rwork = 1. - fwetl(i)
         agcls(i) = rwork * agcls(i)
         agcl3(i) = rwork * agcl3(i)
         agcl4(i) = rwork * agcl4(i)
       if(isimagro .gt. 0) then
         agcc3(i) = rwork * agcc3(i)
         agcc4(i) = rwork * agcc4(i)
       endif
         ancls(i) = rwork * ancls(i)
         ancl3(i) = rwork * ancl3(i)
         ancl4(i) = rwork * ancl4(i)
       if(isimagro .gt. 0) then
         ancc3(i) = rwork * ancc3(i)
         ancc4(i) = rwork * ancc4(i)
       endif
200   continue
      return
end subroutine stomataib
