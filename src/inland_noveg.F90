#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine noveg(kpti, kptj)
! ---------------------------------------------------------------------
!
! if no veg surfaces exist, set prog vars to nominal values
!
! (sensible fluxes fsen[u,s,l], latent fluxes fvap[u,s,l]*, 
! temperature t[u,s,l], and intercepted liquid, snow amounts 
! wliq[u,s,l], wsno[u,s,l] have been calculated for a unit 
! leaf/stem surface, whether or not one exists.)
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsno
      use inland_comsoi
      use inland_comveg

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti ! index of 1st point of little vector in big lpt vector
      integer kptj ! index of last point of little vector

! local variables
      integer i
      real*8 tav, &  ! average temp for soil and snow 
             x,   &  ! total lai + sai
             y       ! fraction of lower canopy not snow covered 

      do 100 i = kpti, kptj
         tav = (1.-fi(i))*tg(i) + fi(i)*ti(i)
         if (lai(i,2).eq.0. .or. fu(i).eq.0.) then
            tu(i) = tav
            wliqu(i) = 0.
            wsnou(i) = 0.
         endif
         if (sai(i,2).eq.0. .or. fu(i).eq.0.) then
            ts(i) = tav
            wliqs(i) = 0.
            wsnos(i) = 0.
         endif 
         x = 2.0 * (lai(i,1) + sai(i,1))
         y = fl(i)*(1.-fi(i))
         if (x .eq.0. .or. y.eq.0.) then
            tl(i) = tav 
            wliql(i) = 0.
            wsnol(i) = 0.
         endif
100   continue
      return
end subroutine noveg
