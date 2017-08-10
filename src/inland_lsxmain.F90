#include "inland_config.h"
! ---------------------------------------------------------------
subroutine lsxmain(loopi, kpti, kptj)
!---------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_comsoi
      use inland_com1d
      use inland_control, only: isinfilt

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer loopi           ! index of little vector in big vector
      integer kpti            ! index of 1st point of little vector 
!                             ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer ib, &           ! waveband number (1=visible, 2= near-IR)
              i               ! loop indice

! Kaiyuan Li for Green-Ampt infiltration model
      if (isinfilt.eq.1)  call initsw(kpti,kptj)

! set physical soil quantities
      call setsoi(kpti, kptj)

! calculate areal fractions wetted by intercepted h2o
      call fwetcal(kpti, kptj)

! set up for solar calculations
      call solset(loopi, kpti, kptj)

! solar calculations for each waveband
      do 100 ib = 1, nband

! solsur sets surface albedos for soil and snow
! solalb performs the albedo calculations
! solarf uses the unit-incident-flux results from solalb
! to obtain absorbed fluxes sol[u,s,l,g,i] and 
! incident pars sunp[u,l]
         call solsur (ib, loopi)
         call solalb (ib, loopi)
         call solarf (ib, loopi)
100   continue

! calculate ir fluxes
      call irrad(kpti, kptj)

! step intercepted h2o
      call cascade(kpti, kptj)

! re-calculate wetted fractions, changed by cascade
      call fwetcal(kpti, kptj)

! step vegetation canopy temperatures implicitly
! and calculate sensible heat and moisture fluxes
      call canopy(kpti, kptj)

! step intercepted h2o due to evaporation
      call cascad2(kpti, kptj)

! arbitrarily set veg temps & intercepted h2o for no-veg locations
      call noveg(kpti, kptj)

! set net surface heat fluxes for soil and snow models
      do 110 i = kpti, kptj 
         heatg(i) = solg(i) + firg(i) - fseng(i) - hvasug(i)*fvapg(i)
         heati(i) = soli(i) + firi(i) - fseni(i) - hvasui(i)*fvapi(i)
 
110   continue

! step snow model
      call snow(kpti, kptj) 

! step soil model
      call soilctl(kpti, kptj)

      return
end subroutine lsxmain
