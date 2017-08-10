#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine vegmap
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comveg
      use inland_comcrop, only:isimagro

      implicit none
!-----------------------------------------------------------------------
!
! local variables
!
      integer i, j

      real*8 treelai, shrublai, grasslai, ratio, &
             maxlai,    & ! maximum lai
             totlai,    & ! total ecosystem lai
             cropbio
!
      integer domtree    ! dominant tree
!-----------------------------------------------------------------------
!
! classify vegetation cover into standard inland vegetation classes 
!
! ---------------------------------------------------
!  1: tropical evergreen forest / woodland
!  2: tropical deciduous forest / woodland
!  3: temperate evergreen broadleaf forest / woodland
!  4: temperate evergreen conifer forest / woodland
!  5: temperate deciduous forest / woodland
!  6: boreal evergreen forest / woodland
!  7: boreal deciduous forest / woodland
!  8: mixed forest / woodland
!  9: savanna
! 10: grassland / steppe 
! 11: dense shrubland
! 12: open shrubland
! 13: tundra
! 14: desert 
! 15: polar desert / rock / ice
! ---------------------------------------------------
!
! begin global grid
!
      do 100 i = lbeg, lend
!
! determine tree, shrub, grass and total lai
!
         treelai  = totlaiu(i) 
         shrublai = plai(i,9)  + plai(i,10)
         grasslai = plai(i,11) + plai(i,12)
    if(isimagro .gt. 0) then
!
! crop biomass -- as used as an overriding condition for
! determining a vegetation class
!
        cropbio  = 0.
        do 105 j = scpft, ecpft
         cropbio   = cropbio + biomass(i,j)
105    continue
    endif
! 
!
         totlai = max (dble(0.01), totlail(i) + totlaiu(i))
!
! determine dominant tree type by lai dominance
!
         domtree = 0
         maxlai = 0.0
!
         do 110 j = 1, 8
            if (plai(i,j).gt.maxlai) then
               domtree = j
               maxlai = plai(i,j)
            endif
 110    continue
!
! assign initial vegetation type
!
        vegtype0(i) = -999.99
!
! dominant type:  tropical broadleaf evergreen tree
!
      if (domtree.eq.1) then
         if (treelai.gt.2.5) vegtype0(i) =  1.0  ! tropical evergreen forest / woodland
         if (treelai.le.2.5) vegtype0(i) =  9.0  ! savanna
         if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
      endif
!
! dominant type:  tropical broadleaf drought-deciduous tree
!
      if (domtree.eq.2) then
         if (treelai.gt.2.5)         vegtype0(i) =  2.0  ! tropical deciduous forest / woodland
         if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
         if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
      endif
!
! dominant type:  warm-temperate broadleaf evergreen tree
!
      if (domtree.eq.3) then
         if (treelai.gt.2.5)         vegtype0(i) =  3.0  ! temperate evergreen broadleaf forest / woodland
         if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
         if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
       endif
!
! dominant type:  temperate conifer evergreen tree
!
      if (domtree.eq.4) then
         if (treelai.gt.1.5)         vegtype0(i) =  4.0  ! temperate evergreen conifer forest / woodland
         if (treelai.le.1.5)         vegtype0(i) =  9.0  ! savanna
         if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
      endif
!
! dominant type:  temperate broadleaf deciduous tree
!
      if (domtree.eq.5) then
         if (treelai.gt.1.5)         vegtype0(i) =  5.0  ! temperate deciduous forest / woodland
         if (treelai.le.1.5)         vegtype0(i) =  9.0  ! savanna
         if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
      endif
!
! dominant type:  boreal conifer evergreen tree
!
      if (domtree.eq.6)             vegtype0(i) =  6.0  ! boreal evergreen forest / woodland
!
!       if (domtree.eq.6) then
!         if (treelai.gt.1.0)         vegtype0(i) =  6.0  ! boreal evergreen forest / woodland
!         if (treelai.le.1.0) then
!           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
!           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
!         endif
!       endif
!
! dominant type:  boreal broadleaf cold-deciduous tree
!
      if (domtree.eq.7) vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
!
!       if (domtree.eq.7) then
!         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
!         if (treelai.le.1.0) then
!           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
!           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
!         endif
!       endif
!
! dominant type:  boreal conifer cold-deciduous tree
!
      if (domtree.eq.8) vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
!
!       if (domtree.eq.8) then
!         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
!         if (treelai.le.1.0) then
!           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
!           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
!         endif
!       endif
!
! temperate/boreal forest mixtures
!
      if ((domtree.ge.4).and.(domtree.le.8)) then
         ratio = (plai(i,5) + plai(i,7) + plai(i,8)) / &
              (plai(i,4) + plai(i,5) + plai(i,6) + &
              plai(i,7) + plai(i,8))
         if (treelai.gt.1.0) then
            if ((ratio.gt.0.45).and.(ratio.lt.0.55)) vegtype0(i) = 8.
         endif
         if ((domtree.le.5).and.(treelai.le.1.0)) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
      endif
!
! no tree is dominant
!
      if (domtree.eq.0) then
         if (treelai.gt.1.0)         vegtype0(i) =  9.0  ! savanna
         if (treelai.le.1.0) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
      endif
!
! overriding vegtation classifications
!
      if (totlai.lt.1.0)            vegtype0(i) = 12.0  ! open shrubland
      if (totlai.le.0.4)            vegtype0(i) = 14.0  ! desert
!
! overriding climatic rules
!
      if (gdd5(i).lt.350.0) then
         if (totlai.ge.0.4)          vegtype0(i) = 13.0  ! tundra
         if (totlai.lt.0.4)          vegtype0(i) = 15.0  ! polar desert
      endif
!
        if (gdd0(i).lt.100.0)         vegtype0(i) = 15.0  ! polar desert
!
        if (cropbio .gt. 0.0)         vegtype0(i) = 16.0  ! croplands
!
 100  continue
!
! return to the main program
!
      return
end subroutine vegmap
