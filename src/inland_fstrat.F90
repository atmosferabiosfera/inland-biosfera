#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine fstrat(tb,tt,ttfac,qb,qt,zb,zt,albm,albh,alt,u,rich,stram,strah, &
                  iter, kpti, kptj)
! ---------------------------------------------------------------------
!
! computes mixing-length stratification correction factors
! for momentum and heat/vapor, for current 1d strip, using
! parameterizations in louis (1979),blm,17,187. first computes
! richardson numbers. sets an upper limit to richardson numbers
! so lower-veg winds don't become vanishingly small in very
! stable conditions (cf, carson and richards,1978,blm,14,68)
!
! system (i) is as in louis(1979). system (vi) is improved as
! described in louis(1982), ecmwf workshop on planetary boundary
! layer parameterizations,november 1981,59-79 (qc880.4 b65w619)
!-----------------------------------------------------------------------
      use inland_parameters

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer iter            ! iteration number
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

      real*8 tb(lbeg:lend),     &  ! bottom temperature (supplied)
             tt(lbeg:lend),     &  ! top temperature (supplied)
             ttfac(lbeg:lend),  &  ! pot. temp factor for ttop (relative to bottom,supplied)
             qb(lbeg:lend),     &  ! bottom specific humidity (supplied)
             qt(lbeg:lend),     &  ! top specific humidity (supplied)
             zb(lbeg:lend),     &  ! height of bottom (supplied)
             zt(lbeg:lend),     &  ! height of top (supplied)
             albm(lbeg:lend),   &  ! log (bottom roughness length) for momentum (supplied)
             albh(lbeg:lend),   &  ! log (bottom roughness length) for heat/h2o (supplied)
             alt(lbeg:lend),    &  ! log (z at top) (supplied)
             u(lbeg:lend),      &  ! wind speed at top (supplied)
             rich(lbeg:lend),   &  ! richardson number (returned)
             stram(lbeg:lend),  &  ! stratification factor for momentum (returned)
             strah(lbeg:lend),  &  ! stratification factor for heat/vap (returned)
             stramx(kpti:kptj), &
             strahx(kpti:kptj)  

! local variables
      integer indp(kptj-kpti+1), indq(kptj-kpti+1)
      integer i, j, np, nq
      real*8 zht, zhb, xm, xh, rwork, ym, yh, z, w
! ---------------------------------------------------------------------
      np = 0
      nq = 0

! do for all points
      do 100 i = kpti, kptj
! calculate richardson numbers
         zht = tt(i)*ttfac(i)*(1.+.622*qt(i))
         zhb = tb(i)*         (1.+.622*qb(i))

         rich(i) = grav * max (zt(i)-zb(i), dble(0.)) * (zht-zhb) / (0.5*(zht+zhb) * u(i)**2)

! bound richardson number between -2.0 (unstable) to 1.0 (stable)
         rich(i) = max (dble(-2.0), min (rich(i), dble(1.0)))
100   continue

! set up indices for points with negative or positive ri
      do 110 i = kpti, kptj
         if (rich(i).le.0.) then
            np = np + 1
            indp(np) = i
         else
            nq = nq + 1
            indq(nq) = i
         endif
110   continue

! calculate momentum and heat/vapor factors for negative ri
      if (np.gt.0) then
         do 200 j = 1, np
            i = indp(j)
            xm = max (alt(i)-albm(i), dble(.5))
            xh = max (alt(i)-albh(i), dble(.5))
            rwork = sqrt(-rich(i))
            ym = (vonk/xm)**2 * exp (0.5*xm) * rwork
            yh = (vonk/xh)**2 * exp (0.5*xh) * rwork

! system (vi)
            stramx(i) = 1.0 - 2*5*rich(i) / (1.0 + 75*ym)
            strahx(i) = 1.0 - 3*5*rich(i) / (1.0 + 75*yh)
200      continue
      endif

! calculate momentum and heat/vapor factors for positive ri
      if (nq.gt.0) then
         do 300 j=1,nq
            i = indq(j)

! system (vi)
            z = sqrt(1.0 + 5 * rich(i))
            stramx(i) = 1.0 / (1.0 + 2*5*rich(i) / z)
            strahx(i) = 1.0 / (1.0 + 3*5*rich(i) * z)
300      continue
      endif

! except for the first iteration, weight results with the
! previous iteration's values. this improves convergence by
! avoiding flip-flop between stable/unstable stratif, eg,
! with cold upper air and the lower surface being heated by
! solar radiation
      if (iter.eq.1) then
         do 400 i = kpti, kptj
            stram(i) = stramx(i)
            strah(i) = strahx(i)
400      continue
      else
         w = 0.5
         do 410 i = kpti, kptj
            stram(i) = w * stramx(i) + (1.0 - w) * stram(i)
            strah(i) = w * strahx(i) + (1.0 - w) * strah(i)
 410     continue
      endif
      return
end subroutine fstrat
