#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine solalb (ib, loopi)
! ---------------------------------------------------------------------
! calculates effective albedos of the surface system,
! separately for unit incoming direct and diffuse flux -- the 
! incoming direct zenith angles are supplied in comatm array 
! coszen, and the effective albedos are returned in comatm
! arrays asurd, asuri -- also detailed absorbed and reflected flux
! info is stored in com1d arrays, for later use by solarf
!
! the procedure is first to calculate the grass+soil albedos,
! then the tree + (grass+soil+snow) albedos. the labels
! (a) to (d) correspond to those in the description doc
!
! ib = waveband number
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d

      implicit none
! ---------------------------------------------------------------------
! input variable
      integer ib
      integer loopi           ! index of little vector in big vector

! local variable
      integer j, &     ! loop indice on number of points with >0 coszen
              nsolz, & ! number of points with >0 coszen
              i

! do nothing if all points in current strip have coszen le 0
      if (nsol(loopi).eq.0) return
      nsolz = nsol(loopi)

! (a) obtain albedos, etc, for two-stream lower veg + soil
!     system, for direct and diffuse incoming unit flux
      do 100 j = 1, nsolz
         i = indsol(loopi,j)
         asurd(i,ib) = albsod(i)
         asuri(i,ib) = albsoi(i)
100   continue

      call twostrib (ablod, abloi,  relod, reloi,  flodd,  dummy, &
                    flodi, floii,  asurd,  asuri,    1,   coszen, ib, &
                    loopi)

! (b) areally average surface albedos (lower veg, soil, snow)
      do 200 j = 1, nsolz
         i = indsol(loopi,j)
         asurd(i,ib) = fl(i)*(1.-fi(i))*relod(i) &
                       + (1.-fl(i))*(1.-fi(i))*albsod(i) &
                       + fi(i)*albsnd(i)    
         asuri(i,ib) = fl(i)*(1.-fi(i))*reloi(i) &
                       + (1.-fl(i))*(1.-fi(i))*albsoi(i) &
                       + fi(i)*albsni(i)
200   continue

! (c) obtain albedos, etc, for two-stream upper veg + surface
!     system, for direct and diffuse incoming unit flux
      call twostrib (abupd, abupi,  reupd, reupi,  fupdd,  dummy, &
                    fupdi, fupii,  asurd,  asuri,    2,   coszen, ib, &
                    loopi)

! (d) calculate average overall albedos 
      do 300 j = 1, nsolz
         i = indsol(loopi,j)
         asurd(i,ib) = fu(i) * reupd(i) + (1.-fu(i)) * asurd(i,ib)
         asuri(i,ib) = fu(i) * reupi(i) + (1.-fu(i)) * asuri(i,ib)
300   continue
      return
end subroutine solalb 
