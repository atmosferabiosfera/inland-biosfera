#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine vadapt (hcur, tcur, htop, indp, np, nlay, kpti, kptj) 
! ---------------------------------------------------------------------
! re-adapt snow layer thicknesses, so top thickness
! equals hsnotop and other thicknesses are equal
!
! also adjusts profile of tracer field tcur so its vertical
! integral is conserved (eg, temperature)
!
! hcur = layer thicknesses (supplied and returned)
! tcur = tracer field (supplied and returned)
! htop = prescribed top layer thickness (supplied)
! indp = index of snow pts in current strip (supplied)
! np   = number of snow pts in current strip (supplied)
      use inland_parameters

      implicit none
!-----------------------------------------------------------------------
! input-output variables variables
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

      integer np, nlay
      real*8 htop

      integer indp(np)
      real*8 hcur(lbeg:lend,nlay),       tcur(lbeg:lend,nlay)         

! local variables
      integer k, i, j, ko
      integer npt             ! number of points in little vector
      real*8 rwork, dz

      real*8 hnew(kpti:kptj,nsnolay),   &  ! new layer thickness
             tnew(kpti:kptj,nsnolay),   &  ! new temperatures of layers 
             zold(kpti:kptj,nsnolay+1), &  ! distances from surface to old layer boundaries
             ht(kpti:kptj),             &  ! storing variable for zold 
             h1(kpti:kptj),             &  ! to compute new layer thickness
             za(kpti:kptj),             &  !
             zb(kpti:kptj),             &
             zheat(kpti:kptj)

!-----------------------------------------------------------------------
! if no snow or seaice points in current 1d strip, return. note
! that the index is not used below (for cray vec and efficiency)
! except in the final loop setting the returned values
      if (np.eq.0) return
         npt = kptj - kpti + 1

! set distances zold from surface to old layer boundaries
      do k = 1, nsnolay+1
         do i = kpti, kptj
            zold(i,k) = 0.0
         end do
      end do

      do 300 k=1,nlay
         do 302 i= kpti, kptj
            zold(i,k+1) = zold(i,k) + hcur(i,k)
302      continue
300   continue

! set new layer thicknesses hnew (tot thickness is unchanged).
! if total thickness is less than nlay*htop (which should be
! le hsnomin), make all new layers equal including
! top one, so other layers aren't so thin. use epsilon to 
! handle zero (snow) points
      do i = kpti, kptj
         ht(i) = zold(i,nlay+1)
      end do

      rwork = nlay*htop
      do 304 i= kpti, kptj
         if (ht(i).ge.rwork) then
            h1(i) = (ht(i)-htop)/(nlay-1)
         else
            h1(i) = max (ht(i)/nlay, epsilon)
         endif
304   continue

      do 306 k=1,nlay
         do 308 i= kpti, kptj
            hnew(i,k) = h1(i)
308      continue
306   continue

      rwork = nlay*htop
      do 310 i= kpti, kptj
         if (ht(i).ge.rwork) hnew(i,1) = htop
310   continue

! integrate old temperature profile (loop 410) over each
! new layer (loop 400), to get new field tnew
      zb(:)=0.0
      do 400 k=1,nlay
         do 402 i= kpti, kptj
            za(i) = zb(i)
            zb(i) = za(i) + hnew(i,k)
402      continue
         zheat(:)=0.0
         do 410 ko=1,nlay
            do 412 i= kpti, kptj
               if (za(i).lt.zold(i,ko+1) .and. zb(i).gt.zold(i,ko)) then
                  dz = min(zold(i,ko+1),zb(i)) - max(zold(i,ko),za(i))
                  zheat(i) = zheat(i) + tcur(i,ko)*dz
               endif
412         continue
410      continue

         do 420 i= kpti, kptj
            tnew(i,k) = zheat(i) / hnew(i,k)
420      continue
400   continue

! use index for final copy to seaice or snow arrays, to avoid
! changing soil values (when called for seaice) and to avoid
! changing nominal snow values for no-snow points (when called
! for snow)
      do 500 k=1,nlay
         do 502 j=1,np
            i = indp(j)
            hcur(i,k) = hnew(i,k)
            tcur(i,k) = tnew(i,k)
502      continue
500   continue
      return
end subroutine vadapt
