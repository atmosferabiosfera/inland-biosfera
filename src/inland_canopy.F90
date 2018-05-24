#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine canopy(kpti, kptj)
! ---------------------------------------------------------------------
!
! calculates sensible heat and moisture flux coefficients,
! and steps canopy temperatures through one timestep
!
! atmospheric conditions at za are supplied in comatm
! arrays ta, qa and psurf
!
! downward sensible heat and moisture fluxes at za
! are returned in com1d arrays fsena, fvapa
!
! sensible heat and moisture fluxes from solid objects to air
! are stored (for other models and budget) in com1d arrays
! fsen[u,s,l,g,i], fvap[u,s,l,g,i]
!
! the procedure is first to compute wind speeds and aerodynamic
! transfer coefficients in turcof, then call turvap to solve an
! implicit linear system for temperatures and specific
! humidities and the corresponding fluxes - this is iterated
! niter times for non-linearities due to stratification,
! implicit/explicit (h2o phase), dew, vpd and max soil
! moisture uptake - t12 and q12 are changed each iteration,
! and tu, ts, tl, tg, ti can be adjusted too
!---------------------------------------------------------------
! common blocks
      use inland_com1d
      use inland_parameters
      use inland_comveg
      use inland_comatm

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector in big lpt vector
      integer kptj            ! index of last point of little vector

!local variable
      integer npt             ! number of points in little vector
      integer niter, iter, i

      real*8 cdmaxa, cdmaxb, ctau

      real*8 xu(lbeg:lend),xs(lbeg:lend),xl(lbeg:lend),chux(lbeg:lend), & 
         chsx(lbeg:lend),chlx(lbeg:lend),chgx(lbeg:lend),wlgx(lbeg:lend), &
         wigx(lbeg:lend),cog(lbeg:lend),coi(lbeg:lend),zirg(lbeg:lend), &
         ziri(lbeg:lend),wu(lbeg:lend),ws(lbeg:lend),wl(lbeg:lend), & 
         wg(lbeg:lend),wi(lbeg:lend),tuold(lbeg:lend),tsold(lbeg:lend), &
         tlold(lbeg:lend),tgold(lbeg:lend),tiold(lbeg:lend),qgfac0(lbeg:lend)

! initialize aerodynamic quantities
      call canini(kpti, kptj)

! estimate soil moisture stress parameters
      call drystress(kpti, kptj)

! iterate the whole canopy physics solution niter times:
      niter = 3

      do 100 iter = 1, niter

! calculate wind speeds and aerodynamic transfer coeffs
         call turcof(iter, kpti, kptj)

! calculate canopy photosynthesis rates and conductance
         call stomataib(kpti, kptj)

! solve implicit system of heat and water balance equations
         call turvap (iter, niter, kpti, kptj, xu, xs, xl, chux, chsx, chlx, &
                      chgx, wlgx, wigx, cog, coi, zirg, ziri, wu, ws, wl, wg, &
                      wi, tuold, tsold, tlold, tgold, tiold, qgfac0)
100   continue

      cdmaxa = 300./(2.*dtime)
      cdmaxb = 1e20
      do i = kpti, kptj
         ctau = ua(i) * (vonk / (aloga(i) - alogu(i)))**2 * stramu(i)
         ctau = min (cdmaxa, ctau / (1. + ctau/cdmaxb))
         taux(i) = rhoa(i) * ctau * ux(i)
         tauy(i) = rhoa(i) * ctau * uy(i)
      end do

! Calculate 2-m surface air temperature (diagnostic, for history)
! Arguments are 1-dimensional --> can be passed only for the kpti-->kptj
! and tscreen doesn't use any array defined over all land points.
      npt = kptj - kpti + 1

      call tscreen  (ts2(kpti), dble(2.), za(kpti), z1(kpti),z12(kpti),z2(kpti), &
                     z3(kpti), z34(kpti), z4(kpti), dispu(kpti),displ(kpti), &
                     ta(kpti),t12(kpti),t34(kpti), npt)

! Also calculate 2-m specific humidity

      call tscreen  (qs2(kpti), dble(2.), za(kpti), z1(kpti), z12(kpti), z2(kpti), &
                     z3(kpti),z34(kpti),z4(kpti), dispu(kpti),displ(kpti), qa(kpti), &
                     q12(kpti),q34(kpti), npt)
      return
end subroutine canopy
