#include "inland_config.h"
module inland_com1d

      implicit none
      public
      save
!
! -----
! com1d
! -----
!
! canopy scaling terms (Ramankutty, 1998 in progress)
!     terml(npoi,7)      ! term needed in lower canopy scaling
!     termu(npoi,7)      ! term needed in upper canopy scaling
      real*8, dimension(:,:), allocatable :: terml, termu

!     scalcoefl(npoi,4)  ! term needed in lower canopy scaling
!     scalcoefu(npoi,4)  ! term needed in upper canopy scaling
      real*8, dimension(:,:), allocatable :: scalcoefl, scalcoefu

! variables for internal fluxes: passed to soil and snow models

!     firb (npoi)        ! net upward ir radiation at reference atmospheric level za (W m-2)
!     firg (npoi)        ! ir radiation absorbed by soil/ice (W m-2)
!     firi (npoi)        ! ir radiation absorbed by snow (W m-2)
!     firl (npoi)        ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
!     firs (npoi)        ! ir radiation absorbed by upper canopy stems (W m-2)
!     firu (npoi)        ! ir raditaion absorbed by upper canopy leaves (W m-2)
!     fsena (npoi)       ! downward sensible heat flux between za & z12 at za (W m-2)
!     fseng (npoi)       ! upward sensible heat flux between soil surface & air at z34 (W m-2)
!     fseni(npoi)        ! upward sensible heat flux between snow surface & air at z34 (W m-2)
!     fsenl (npoi)       ! sensible heat flux from lower canopy to air (W m-2)
!     fsens (npoi)       ! sensible heat flux from upper canopy stems to air (W m-2)
!     fsenu (npoi)       ! sensible heat flux from upper canopy leaves to air (W m-2)
!     fvapa (npoi)       ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
!     fvapg (npoi)       ! h2o vapor flux (evaporation) between soil & air at z34 (kg m-2 s-1/bare ground fraction)
!     fvapi (npoi)       ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
!     fvaplt (npoi)      ! h2o vapor flux (transpiration) between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
!     fvaplw (npoi)      ! h2o vapor flux (evaporation from wet surface) between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
!     fvaps (npoi)       ! h2o vapor flux (evaporation from wet surface) between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
!     fvaput (npoi)      ! h2o vapor flux (transpiration from dry parts) between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
!     fvapuw (npoi)      ! h2o vapor flux (evaporation from wet parts) between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
!     soli (npoi)        ! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
!     solg (npoi)        ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
!     soll (npoi)        ! solar flux (direct + diffuse) absorbed by lower canopy leaves and stems per unit canopy area (W m-2)
!     sols (npoi)        ! solar flux (direct + diffuse) absorbed by upper canopy stems per unit canopy area (W m-2)
!     solu (npoi)        ! solar flux (direct + diffuse) absorbed by upper canopy leaves per unit canopy area (W m-2)
!     raing(npoi)        ! rainfall rate at soil level (kg m-2 s-1)
!     rainl(npoi)        ! rainfall rate below upper canopy (kg m-2 s-1)
!     rainu(npoi)        ! rainfall rate above upper canopy (kg m-2 s-1)
!     snowg(npoi)        ! snowfall rate at soil level (kg h2o m-2 s-1)
!     snowl(npoi)        ! snowfall rate below upper canopy (kg h2o m-2 s-1)
!     snowu(npoi)        ! snowfall rate above upper canopy (kg h2o m-2 s-1)
!     traing(npoi)       ! rainfall temperature at soil level (K)
!     trainl(npoi)       ! rainfall temperature below upper canopy (K)
!     trainu(npoi)       ! rainfall temperature above upper canopy (K)
!     tsnowg(npoi)       ! snowfall temperature at soil level (K) 
!     tsnowl(npoi)       ! snowfall temperature below upper canopy (K)
!     tsnowu(npoi)       ! snowfall temperature above upper canopy (K)
!     flx(20)            ! sum dtime flux output into daily value
      real*8, dimension(:), allocatable :: firb, firg, firi, firl, firs, firu, fsena, fseng, fseni, &
                                           fsenl, fsens, fsenu, fvapa, fvapg, fvapi, fvaplt, fvaplw, &
                                           fvaps, fvaput, fvapuw, soli, solg, soll, sols, solu, raing, &
                                           rainl, rainu, snowg, snowl, snowu, traing, trainl, trainu, &
                                           tsnowg, tsnowl, tsnowu
      real*8  flx(20)
!
! variables for solar calculations
!
! note that all direct fluxes are per unit horizontal area
! (i.e., already including a factor of cos(zen angle))
!     nsol(numlv)        ! number of points in indsol
      integer, dimension(:), allocatable :: nsol

!     abupd(npoi)        ! fraction of direct  radiation absorbed by upper canopy
!     abupi(npoi)        ! fraction of diffuse radiation absorbed by upper canopy
!     ablod(npoi)        ! fraction of direct  radiation absorbed by lower canopy
!     abloi(npoi)        ! fraction of diffuse radiation absorbed by lower canopy
!     albsnd(npoi)       ! direct  albedo for snow surface (visible or IR)
!     albsni(npoi)       ! diffuse albedo for snow surface (visible or IR)
!     albsod(npoi)       ! direct  albedo for soil surface (visible or IR)
!     albsoi(npoi)       ! diffuse albedo for soil surface (visible or IR)
!     dummy(npoi)        ! placeholder,  always = 0: no direct flux produced for diffuse incident
!     flodd(npoi)        ! downward direct radiation per unit incident direct radiation on lower canopy (W m-2)
!     flodi(npoi)        ! downward diffuse radiation per unit incident direct radiation on lower canopy (W m-2)
!     floii(npoi)        ! downward diffuse radiation per unit incident diffuse radiation on lower canopy
!     fupdd(npoi)        ! downward direct radiation per unit incident direct beam on upper canopy (W m-2)
!     fupdi(npoi)        ! downward diffuse radiation per unit icident direct radiation on upper canopy (W m-2)
!     fupii(npoi)        ! downward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
!     relod(npoi)        ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
!     reloi(npoi)        ! upward diffuse radiation per unit incident diffuse radiation on lower canopy (W m-2)
!     reupd(npoi)        ! upward direct radiation per unit incident direct radiation on upper canopy (W m-2)
!     reupi(npoi)        ! upward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
!     sol2d(npoi)        ! direct downward radiation  out of upper canopy per unit vegetated (upper) area (W m-2)
!     sol2i(npoi)        ! diffuse downward radiation out of upper canopy per unit vegetated (upper) area(W m-2)
!     sol3d(npoi)        ! direct downward radiation  out of upper canopy + gaps per unit grid cell area (W m-2)
!     sol3i(npoi)        ! diffuse downward radiation out of upper canopy + gaps per unit grid cell area (W m-2)
      real*8, dimension(:), allocatable :: abupd, abupi, ablod, abloi, albsnd, albsni, albsod, albsoi, &
                                           dummy, flodd, flodi, floii, fupdd, fupdi, fupii, relod, reloi, &
                                           reupd, reupi, sol2d, sol2i, sol3d, sol3i

!     indsol(numlv,npoi) ! index of current strip for points with positive coszen
      integer, dimension(:,:), allocatable :: indsol

! variables for aerodynamic calculations
!
!     aloga(npoi)        ! log (za - dispu) 
!     alogav(npoi)       ! average of alogi and alogg 
!     alogg(npoi)        ! log of soil roughness
!     alogi(npoi)        ! log of snow roughness
!     alogl(npoi)        ! log (roughness length of lower canopy)
!     alogu(npoi)        ! log (roughness length of upper canopy)
!     alog1(npoi)        ! log (z1 - dispu) 
!     alog2(npoi)        ! log (z2 - displ)
!     alog3(npoi)        ! log (z3 - displ)
!     alog4(npoi)        ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
!     bdl(npoi)          ! aerodynamic coefficient ([(tau/rho)/u**2] for laower canopy (A31/A30 Pollard & Thompson 1995)
!     bdu(npoi)          ! aerodynamic coefficient ([(tau/rho)/u**2] for upper canopy (A31/A30 Pollard & Thompson 1995)
!     cl(npoi)           ! air transfer coefficient (*rhoa) (m s-1 kg m-3) between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
!     cp(npoi)           ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
!     cu(npoi)           ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
!     dil(npoi)          ! inverse of momentum diffusion coefficient within lower canopy (m)
!     displ(npoi)        ! zero-plane displacement height for lower canopy (m)
!     dispu(npoi)        ! zero-plane displacement height for upper canopy (m)
!     diu(npoi)          ! inverse of momentum diffusion coefficient within upper canopy (m)
!     expl(npoi)         ! exphl**2
!     expu(npoi)         ! exphu**2
!     exphl(npoi)        ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
!     exphu(npoi)        ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
!     pfluxl(npoi)       ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
!     pfluxs(npoi)       ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
!     pfluxu(npoi)       ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
!     rhoa(npoi)         ! air density at za (allowing for h2o vapor) (kg m-3)
!     richl(npoi)        ! richardson number for air above upper canopy (z3 to z2)
!     richu(npoi)        ! richardson number for air between upper & lower canopy (z1 to za)
!     sg(npoi)           ! air-soil transfer coefficient
!     si(npoi)           ! air-snow transfer coefficient
!     strahl(npoi)       ! heat/vap correction factor for stratif between upper & lower canopy (z3 to z2) (louis et al.)
!     strahu(npoi)       ! heat/vap correction factor for stratif above upper canopy (z1 to za) (louis et al.)
!     straml(npoi)       ! momentum correction factor for stratif between upper & lower canopy (z3 to z2) (louis et al.)
!     stramu(npoi)       ! momentum correction factor for stratif above upper canopy (z1 to za) (louis et al.)
!     u1(npoi)           ! wind speed at level z1 (m s-1)
!     u12(npoi)          ! wind speed at level z12 (m s-1)
!     u2(npoi)           ! wind speed at level z2 (m s-1)
!     u3(npoi)           ! wind speed at level z3 (m s-1)
!     u34(npoi)          ! wind speed at level z34 (m s-1)
!     u4(npoi)           ! wind speed at level z4 (m s-1)
!     za(npoi)           ! height above the surface of atmospheric forcing (m)
!     z1(npoi)           ! effective top of upper canopy (for momentum) (m)
!     z12(npoi)          ! effective middle of the upper canopy (for momentum) (m)
!     z2(npoi)           ! effective bottom of the upper canopy (for momentum) (m)
!     z3(npoi)           ! effective top of the lower canopy (for momentum) (m)
!     z34(npoi)          ! effective middle of the lower canopy (for momentum) (m)
!     z4(npoi)           ! effective bottom of the lower canopy (for momentum) (m)
!     taux(npoi)         ! eastward wind stress to sfc(below ZA) from atmos
!     tauy(npoi)         ! poleward wind stress to sfc(below ZA) from atmos
!     ts2(npoi)          ! 2-m surface-air temperature
!     qs2(npoi)          ! 2-m specific humidity
!
!     tfaca(npoi)        ! (ps/p) ** (rair/cair) for atmospheric level  (const)
      real*8, dimension(:), allocatable :: aloga, alogav, alogg, alogi, alogl, alogu, alog1, alog2, alog3, alog4, &
                                           bdl, bdu, cl, cp, cu, dil, displ, dispu, diu, expl, expu, exphl, exphu, &
                                           pfluxl, pfluxs, pfluxu, rhoa, richl, richu, sg, si, strahl, strahu, straml, &
                                           stramu, u1, u12, u2, u3, u34, u4, za, z1, z12, z2, z3, z34, z4, taux, tauy, &
                                           ts2, qs2, tfaca

! variables for intercepted water
!
!     fwetl (npoi)       ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
!     fwets (npoi)       ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
!     fwetu (npoi)       ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
!     fwetlx (npoi)      ! fraction of lower canopy leaf and stem area wetted if dew forms
!     fwetsx (npoi)      ! fraction of upper canopy stem area wetted if dew forms
!     fwetux (npoi)      ! fraction of upper canopy leaf area wetted if dew forms
!     rliql (npoi)       ! proportion of fwetl due to liquid
!     rliqs (npoi)       ! proportion of fwets due to liquid
!     rliqu (npoi)       ! proportion of fwetu due to liquid
      real*8, dimension(:), allocatable :: fwetl, fwets, fwetu, fwetlx, fwetsx, fwetux, rliql, rliqs, rliqu
!

end module inland_com1d
