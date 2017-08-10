#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine impexp (wimp, tveg, ch, wliq, wsno, iter, kpti, kptj)
! ---------------------------------------------------------------------
! sets the implicit vs explicit fraction in turvap calcs for
! upper leaves, upper stems or lower veg. this is to account for
! temperatures of freezing/melting intercepted h2o constrained
! at the melt point. if a purely implicit calc is used for such
! a surface, the predicted temperature would be nearly the atmos
! equil temp with little sensible heat input, so the amount of
! freezing or melting is underestimated. however, if a purely
! explicit calc is used with only a small amount of intercepted
! h2o, the heat exchange can melt/freeze all the h2o and cause
! an unrealistic huge change in the veg temp. the algorithm
! below attempts to avoid both pitfalls
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
! ---------------------------------------------------------------------
! input/output variables
      integer iter ! current iteration number (supplied)
      integer kpti ! index of 1st point of little vector in big lpt vector
      integer kptj ! index of last point of little vector

      real*8 wimp(lbeg:lend), & ! implicit/explicit fraction (0 to 1) (returned)
             tveg(lbeg:lend), & ! temp of veg (previous iteration's soln) (supp)
             ch(lbeg:lend),   & ! heat capacity of veg (supplied)
             wliq(lbeg:lend), & ! veg intercepted liquid (supplied)
             wsno(lbeg:lend)    ! veg intercepted snow (supplied)

! local variables
      integer i
      integer npt             ! number of points in little vector

      real*8 h, z, winew

! for first iteration, set wimp to fully implicit, and return
      npt = kptj - kpti + 1
      if (iter.eq.1) then
         wimp(:)=1.0
         return
      endif

! for second and subsequent iterations, estimate wimp based on
! the previous iterations's wimp and its resulting tveg.
!
! calculate h, the "overshoot" heat available to melt any snow
! or freeze any liquid. then the explicit fraction is taken to
! be the ratio of h to the existing h2o's latent heat (ie, 100%
! explicit calculation if not all of the h2o would be melted or
! frozen). so winew, the implicit amount, is 1 - that ratio.
! but since we are using the previous iteration's t* results
! for the next iteration, to ensure convergence we need to damp
! the returned estimate wimp by averaging winew with the 
! previous estimate. this works reasonably well even with a
! small number of iterations (3), since for instance with large
! amounts of h2o so that wimp should be 0., a good amount of 
! h2o is melted or frozen with wimp = .25
      do 100 i = kpti, kptj
         h = ch(i) * (tveg(i) - tmelt)
         z = max (abs(h), epsilon)
         winew = 1.0
         if (h.gt.epsilon)  winew = 1. - min (dble(1.), hfus * wsno(i) / z)
         if (h.lt.-epsilon) winew = 1. - min (dble(1.), hfus * wliq(i) / z)
         wimp(i) = 0.5 * (wimp(i) + winew)
100   continue
      return
end subroutine impexp
