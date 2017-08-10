#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine solarf(ib, loopi)
! ---------------------------------------------------------------------
!
! calculates solar fluxes absorbed by upper and lower stories,
! soil and snow
!
! zenith angles are in comatm array coszen, and must be the same
! as supplied earlier to solalb
!
! solarf uses the results obtained earlier by solalb and 
! stored in com1d arrays. the absorbed fluxes are returned in
! com1d arrays sol[u,s,l,g,i]
!
! the procedure is first to calculate the upper-story absorbed
! fluxes and fluxes below the upper story, then the lower-story
! absorbed fluxes and fluxes below the lower story, then fluxes
! absorbed by the soil and snow
!
! ib = waveband number
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d

      implicit none
! ---------------------------------------------------------------------
!
! input variable
      integer ib           ! waveband number (1= visible, 2= near-IR)
      integer loopi        ! index of little vector in big vector
!
! local variables
      integer j,     & ! loop indice on number of points with >0 coszen
              nsolz, & ! number of points with >0 coszen
              i
      real*8 x, y, xd, xi, &
             xaiu,   & ! total single-sided lai+sai, upper canopy
             xail      ! total single-sided lai+sai, lower canopy
!
! do nothing if all points in current strip have coszen le 0
      if (nsol(loopi).eq.0) return

      nsolz = nsol(loopi)

! (f) calculate fluxes absorbed by upper leaves and stems,
!     and downward fluxes below upper veg, using unit-flux
!     results of solalb(c) (apportion absorbed flux between
!     leaves and stems in proportion to their lai and sai)
      do 600 j=1,nsolz
         i = indsol(loopi,j)
         x = solad(i,ib)*abupd(i) + solai(i,ib)*abupi(i)
         y = lai(i,2) / max (lai(i,2)+sai(i,2), epsilon)
         solu(i) = solu(i) + x * y
         sols(i) = sols(i) + x * (1.-y)
         sol2d(i) = solad(i,ib)*fupdd(i)
         sol2i(i) = solad(i,ib)*fupdi(i) + solai(i,ib)*fupii(i)
600   continue

! (g) areally average fluxes to lower veg, soil, snow
      do 700 j=1,nsolz
         i = indsol(loopi,j)
         sol3d(i) = fu(i)*sol2d(i) + (1.-fu(i))*solad(i,ib)
         sol3i(i) = fu(i)*sol2i(i) + (1.-fu(i))*solai(i,ib)
700   continue

! (h,i) calculate fluxes absorbed by lower veg, snow-free soil
!       and snow, using results of (g) and unit-flux results
!       of solalb(a)
      do 800 j=1,nsolz
         i = indsol(loopi,j)
         soll(i) = soll(i) + sol3d(i)*ablod(i) + sol3i(i)*abloi(i)

         xd = (fl(i)*flodd(i) + 1.-fl(i)) * sol3d(i)
         xi = fl(i)*(sol3d(i)*flodi(i)+sol3i(i)*floii(i))+(1.-fl(i))*sol3i(i)

         solg(i) = solg(i) + (1.-albsod(i))*xd + (1.-albsoi(i))*xi
	     solsoi(i) = solsoi(i) + xd + xi
         soli(i) = soli(i) + (1.-albsnd(i))*sol3d(i) + (1.-albsni(i))*sol3i(i)
800   continue

! estimate absorbed pars at top of canopy, toppar[u,l] and
! some canopy scaling parameters
!
! this neglects complications due to differing values of dead vs 
! live elements, averaged into rhoveg, tauveg in vegdat, and 
! modifications of omega due to intercepted snow in twoset
!
! do only for visible band (ib=1)
      if (ib.eq.1) then
         do 900 j = 1, nsolz
            i = indsol(loopi,j)

! the canopy scaling algorithm assumes that the net photosynthesis
! is proportional to absored par (apar) during the daytime. during night,
! the respiration is scaled using a 10-day running-average daytime canopy
! scaling parameter.
!
! apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
!
! some of the required terms (i.e. term[u,l] are calculated in the subroutine 'twostrib'.
! in the equations below, 
!
!   A = scalcoefu(i,1) = term[u,l](i,1) * ipardir(0)
!   B = scalcoefu(i,2) = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
!   C = scalcoefu(i,3) = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
!   A + B + C = scalcoefu(i,4) = also absorbed par at canopy of canopy by leaves & stems
!
! upper canopy:
!
! total single-sided lai+sai
            xaiu = max (lai(i,2)+sai(i,2), epsilon)

! some terms required for use in canopy scaling:
            scalcoefu(i,1) = termu(i,1) * solad(i,ib)
            scalcoefu(i,2) = termu(i,2) * solad(i,ib) + termu(i,3) * solai(i,ib)
            scalcoefu(i,3) = termu(i,4) * solad(i,ib) + termu(i,5) * solai(i,ib)
            scalcoefu(i,4) = scalcoefu(i,1) + scalcoefu(i,2) + scalcoefu(i,3)

! apar of the "top" leaves of the canopy
            topparu(i) = scalcoefu(i,4) * lai(i,2) / xaiu

! lower canopy:
!
! total single-sided lai+sai
            xail  = max (lai(i,1)+sai(i,1), epsilon)
            scalcoefl(i,1) = terml(i,1) * sol3d(i)
            scalcoefl(i,2) = terml(i,2) * sol3d(i) + terml(i,3) * sol3i(i)
            scalcoefl(i,3) = terml(i,4) * sol3d(i) + terml(i,5) * sol3i(i)
            scalcoefl(i,4) = scalcoefl(i,1) + scalcoefl(i,2) + scalcoefl(i,3)

! apart of the "top" leaves of the canopy
            topparl(i) = scalcoefl(i,4) * lai(i,1) / xail
900      continue
      endif
      return
end subroutine solarf
