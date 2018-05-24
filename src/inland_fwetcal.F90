#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine fwetcal(kpti, kptj)
! ---------------------------------------------------------------------
!
!   Calculates fwet[u,s,l], the fractional areas wetted by
! intercepted h2o (liquid and snow combined) -  the maximum value
! fmax (<1) allows some transpiration even in soaked conditions
!
!   Use a linear relation between fwet* and wliq*,wsno* (at least
! for small values), so that the implied "thickness" is constant
! (equal to wliq*max, wsno*max as below) and the typical amount
! evaporated in one timestep in steph2o will not make wliq*,wsno*
! negative and thus cause inland_a spurious unrecoverable h2o loss
! (the max(w*max,.01) below numericaly allows w*max = 0 without
! blowup.) in fact evaporation in one timestep *does* sometimes
! exceed wliq*max (currently 1 kg/m2), so there is an additional
! safeguard in turvap that limits the wetted-area aerodynamic
! coefficients suw,ssw,slw -- if that too fails, there is an
! ad-ho! adjustment in steph2o2 to reset negative wliq*,wsno*
! amounts to zero by taking some water vapor from the atmosphere.
!
!   Also sets rliq[u,s,l], the proportion of fwet[u,s,l] due to
! liquid alone. fwet,rliq are used in turvap, rliq in steph2o.
! (so rliq*fwet, (1-rliq)*fwet are the fractional areas wetted
! by liquid and snow individually.) if fwet is 0, choose rliq
! = 1 if t[u,s,l] ge tmelt or 0 otherwize, for use inland_by turvap and
! steph2o in case of initial dew formation on dry surface.
! ---------------------------------------------------------------------
      use inland_com1d
      use inland_parameters
      use inland_control
      use inland_comveg
      use inland_comcrop, only:isimagro

      implicit none
!-----------------------------------------------------------------------
!
! Input arguments
      integer :: kpti            ! index of 1st point of little vector
                                 ! in big lpt vector
      integer :: kptj            ! index of last point of little vector

! local variables
      integer :: i      ! loop indice
      real*8 :: fmax, & ! maximum water cover on two-sided leaf
                xliq, & ! fraction of wetted leaf (liquid only)
                xtot    ! fraction of wetted leaf (liquid and snow)

! maximum water cover on two-sided leaf
!
!      parameter (fmax = 0.08)
      if(isimagro .gt. 0)then
         fmax = 0.08
      else
         fmax = 0.25
      endif
!
! upper leaves
      do 100 i = kpti, kptj
         xliq = wliqu(i) / max (wliqumax,dble(0.01))
         xtot = xliq + wsnou(i) / max (wsnoumax,dble(0.01))
         fwetu(i) = min (fmax, xtot)
         rliqu(i) = xliq / max (xtot, epsilon)

         if (fwetu(i) .eq. 0.0) then
            rliqu(i) = 1.0
            if (tu(i) .lt. tmelt) rliqu(i) = 0.0
         endif
100   continue

! upper stems
      do 200 i = kpti, kptj
         xliq = wliqs(i) / max(wliqsmax,dble(0.01))
         xtot = xliq + wsnos(i) / max(wsnosmax,dble(0.01))
         fwets(i) = min (fmax, xtot)
         rliqs(i) = xliq / max (xtot, epsilon)

         if (fwets(i) .eq. 0.0) then
            rliqs(i) = 1.0
            if (ts(i) .lt. tmelt) rliqs(i) = 0.0
         endif
200   continue

! lower veg
      do 300 i = kpti, kptj
         xliq = wliql(i) / max(wliqlmax,dble(0.01))
         xtot = xliq + wsnol(i) / max(wsnolmax,dble(0.01))
         fwetl(i) = min (fmax, xtot)

         rliql(i) = xliq / max (xtot, epsilon)
         if (fwetl(i) .eq. 0.0) then
            rliql(i) = 1.0
            if (tl(i) .lt. tmelt) rliql(i) = 0.0
         endif
300   continue
      return
end subroutine fwetcal
