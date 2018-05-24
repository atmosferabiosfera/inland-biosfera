#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine solset(loopi, kpti, kptj)
! ---------------------------------------------------------------------
!  zeros albedos and internal absorbed solar fluxes, and sets
!  index for other solar routines. the index indsol, with number
!  of points nsol, points to current 1d strip arrays whose coszen 
!  values are gt 0 (indsol, nsol are in com1d)
! -----------------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_comsno
      use inland_comveg
      use inland_com1d

      implicit none
! ---------------------------------------------------------------------
! input variables
      integer loopi           ! index of little vector in big vector
      integer kpti            ! index of 1st point of little vector 
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer npt      ! number of points in little vector
      integer j, &      ! loop indice on number of points with >0 coszen
              nsolz, &   ! number of points with >0 coszen
              i

      npt = kptj - kpti + 1

      do j = 1, nband
         do i = kpti, kptj
            asurd(i,j) = 0.0
            asuri(i,j) = 0.0
         end do
      end do

! zeros absorbed solar fluxes sol[u,s,l,g,i]1 since only points
! with +ve coszen will be set in solarf, and since
! sol[u,l,s,g,i]1 are summed over wavebands in solarf
!
! similarly zero par-related arrays set in solarf for turvap
!      solu(:)=0.0
!      sols(:)=0.0
!      soll(:)=0.0
!      solg(:)=0.0
!      soli(:)=0.0
!      topparu(:)=0.0
!      topparl(:)=0.0
      solu(kpti:kptj)=0.0
      sols(kpti:kptj)=0.0
      soll(kpti:kptj)=0.0
      solg(kpti:kptj)=0.0
      soli(kpti:kptj)=0.0
      topparu(kpti:kptj)=0.0
      topparl(kpti:kptj)=0.0

! set canopy scaling coefficients for night-time conditions
      do j = 1, 4
         do i = kpti, kptj
            scalcoefl(i,j) = 0.0
            scalcoefu(i,j) = 0.0
         end do
      end do

! set index of points with positive coszen
      nsolz = 0

      do 300 i = kpti, kptj
         if (coszen(i).gt.0.) then
            nsolz = nsolz + 1
            indsol(loopi, nsolz) = i
         endif
300   continue

      nsol(loopi) = nsolz
      return
end subroutine solset
