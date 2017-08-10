#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine snow(kpti, kptj)
! ---------------------------------------------------------------------
! steps snow model through one timestep
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d
      use inland_comcrop, only:isimagro

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local arrays:
!
!     hinit  = initial layer thicknesses when snow first forms
!     indsno = index of points with snow in current 1d strip
!     fiold  = old fi at start of this timestep
!     fhtop  = heat flux into upper snow surface
!     sflo   = heat flux across snow and buried-lower-veg layer bdries
!     zmelt  = liquid mass flux increments to soil, at temperature 
!              tmelt, due to processes occuring during this step
!     zheat  = heat flux to soil, due to processes occuring this step
!     dfi    = change in fi
!     xl     = lower veg density
!     x*,ht  = temporary arrays
      integer i, k, npt, &
              npn                  ! index indsno, npcounter for pts with snow 

      real*8 rwork, rwork2, &
             finew,         &      ! storing variable for fi
             zhh,           &      ! 0.5*hsnomin
             zdh                   ! max height of snow above hsnomin (?) 

      integer indsno(kptj-kpti+1)

! hinit(nsnolay): initial layer thicknesses when snow first forms
! hsnoruf(kpti:kptj): heigth of snow forced to cover lower canopy (?)
! fiold(kpti:kptj): old fi at start of this timestep
! fhtop(kpti:kptj): heat flux into upper snow surface
! sflo(kpti:kptj,nsnolay+2): heat flux across snow and buried-lower-veg 
!                           layer boundaries
! zmelt(kpti:kptj): liquid mass flux increments to soil, at temperature 
!                  tmelt, due to processes occuring during this step
! zheat(kpti:kptj): heat flux to soil, due to processes occuring this step
! dfi(kpti:kptj): change in fi
! xl(kpti:kptj): lower veg density
! xh, xm, ht, x1, x2, x3 (all kpti:kptj): temporary arrays 
      real*8 hinit(nsnolay), hsnoruf(kpti:kptj), fiold(kpti:kptj), &
             fhtop(kpti:kptj), sflo(kpti:kptj,nsnolay+2), zmelt(kpti:kptj), &
             zheat(kpti:kptj), dfi(kpti:kptj), xl(kpti:kptj), xh(kpti:kptj), &
             xm(kpti:kptj), ht(kpti:kptj), x1(kpti:kptj), x2(kpti:kptj), &
             x3(kpti:kptj)

      do 10 i = kpti, kptj
        if(isimagro .eq. 0) then
            hsnoruf(i) =  min (dble(0.70), max (hsnomin+dble(.05), fl(i)*ztop(i,1)))
        else
            hsnoruf(i) = max (0.10, hsnomin+0.01)
        endif
         xl(i) = fl(i) * 2.0 * (lai(i,1) + sai(i,1))
         x1(i) = tlsub(i)
10    continue
      hinit(1) = hsnotop
      do 15 k = 2, nsnolay
         hinit(k) = (hsnomin - hsnotop) / (nsnolay-1)
15    continue
      do 20 i = kpti, kptj
         fiold(i) = fi(i)
20    continue

! zero out arrays
      npt = kptj - kpti + 1

      do k = 1, nsnolay+2
         do i = kpti, kptj
            sflo(i,k) = 0
         end do
      end do
      zmelt(:)=0.0
      zheat(:)=0.0

! set up index indsno, npn for pts with snow - indsno is used
! only by vadapt - elsewhere below, just test on npn > 0
      npn = 0                                        

      do 30 i = kpti, kptj 
         if (fi(i).gt.0.) then
            npn = npn + 1
            indsno(npn) = i
         endif 
30    continue

! set surface heat flux fhtop and increment top layer thickness
! due to rainfall, snowfall and sublimation on existing snow
      if (npn.gt.0) then
         rwork = dtime / rhos
         do 40 i = kpti, kptj
            fhtop(i) = heati(i) + rainl(i)*(                                   &
                          ch2o*(trainl(i)-tmelt)+hfus + cice*(tmelt-tsno(i,1)) &
                       ) + snowl(i)*cice*(tsnowl(i)-tsno(i,1))
            if (fi(i).gt.0.) &
               hsno(i,1)=hsno(i,1)+(rainl(i)+snowl(i)-fvapi(i))*rwork
40       continue
      endif

! step temperatures due to heat conduction, including buried
! lower-veg temperature tlsub
      if (npn.gt.0) then
         call scopya (npt, tlsub(kpti), x1(kpti))
         call snowheat (tlsub, fhtop, sflo, xl, chl, kpti, kptj)
      endif

! put snowfall from 1-fi snow-free area onto side of existing
! snow, or create new snow if current fi = 0. also reset index.
! (assumes total depth of newly created snow = hsnomin.)
! (fi will not become gt 1 here if one timestep's snowfall
! <= hsnomin, but protect against this anyway.)
!
! if no adjacent snowfall or fi = 1, dfi = 0, so no effect
      ht(:)=0.0
      do 190 k=1,nsnolay
         do 192 i=kpti,kptj
            ht(i) = ht(i) + hsno(i,k)
192      continue
190   continue

      do 195 i=kpti,kptj
         if (ht(i).eq.0.) ht(i) = hsnomin
195   continue

      rwork = dtime / rhos
      do 200 i=kpti,kptj
         dfi(i) = (1.-fi(i))*rwork*snowg(i) / ht(i)
         dfi(i) = min (dfi(i), 1.-fi(i))
200   continue

      do 210 k=1,nsnolay
         do 212 i=kpti,kptj
            if (fi(i)+dfi(i).gt.0.) &
               tsno(i,k) = (tsno(i,k)*fi(i) + tsnowg(i)*dfi(i))/(fi(i)+dfi(i))

! set initial thicknesses for newly created snow
            if (fi(i).eq.0. .and. dfi(i).gt.0.) hsno(i,k) = hinit(k)
212      continue
210   continue

      npn = 0
      do 220 i=kpti,kptj
         fi(i) = fi(i) + dfi(i)
         if (fi(i).gt.0.) then 
            npn = npn + 1
            indsno(npn) = i
         endif
220   continue

! melt from any layer (due to implicit heat conduction, any
! layer can exceed tmelt, not just the top layer), and reduce
! thicknesses (even to zero, and give extra heat to soil)
!
! ok to do it for non-snow points, for which xh = xm = 0
      if (npn.gt.0) then
         rwork = 1. / rhos
         do 300 k=1,nsnolay
            do 302 i=kpti,kptj
               xh(i) = rhos*hsno(i,k)*cice * max(tsno(i,k)-tmelt, dble(0.))
               xm(i) = min (rhos*hsno(i,k), xh(i)/hfus)
               hsno(i,k) = hsno(i,k) - xm(i)*rwork
               tsno(i,k) = min (tsno(i,k),tmelt)
               zmelt(i) = zmelt(i) + fi(i)*xm(i)
               zheat(i) = zheat(i) + fi(i)*(xh(i)-hfus*xm(i))
302         continue
300      continue

! adjust fi and thicknesses for coverage-vs-volume relation
! ie, total thickness = hsnomin for fi < fimax, and fi <= fimax.
! (ok to do it for no-snow points, for which ht=fi=finew=0.)
         ht(:)=0.0
         do 400 k=1,nsnolay
            do 402 i=kpti,kptj
               ht(i) = ht(i) + hsno(i,k)
402         continue
400      continue

! linear variation  for 0 < fi < 1
         zhh = 0.5*hsnomin
         do 404 i=kpti,kptj
            zdh = hsnoruf(i)-hsnomin
            finew = ( -zhh + sqrt(zhh**2 + zdh*fi(i)*ht(i)) ) / zdh
            finew = max (dble(0.), min (fimax, finew))
            x1(i) =  fi(i) / max (finew, epsilon)
            fi(i) =  finew
404      continue
         do 406 k=1,nsnolay
            do 408 i=kpti,kptj
               hsno(i,k) = hsno(i,k) * x1(i)
408         continue
406      continue
      endif

! re-adapt snow thickness profile, so top thickness = hsnotop
! and other thicknesses are equal
!
! adjust temperature to conserve sensible heat
      call vadapt(hsno,tsno,hsnotop,indsno,npn,nsnolay,kpti,kptj)

! if fi is below fimin, melt all snow and adjust soil fluxes
      if (npn.gt.0) then
         call scopya (npt, fi(kpti), x1(kpti))
         do 500 k=1,nsnolay
            do 502 i=kpti,kptj
               if (x1(i).lt.fimin) then
                  xm(i) = x1(i) * rhos * hsno(i,k)
                  zmelt(i) = zmelt(i) + xm(i)
                  zheat(i) = zheat(i) - xm(i)*(cice*(tmelt-tsno(i,k))+hfus)
                  hsno(i,k) = 0.
                  tsno(i,k) = tmelt
                  fi(i) = 0.
               endif
502         continue
500      continue
      endif

! adjust buried lower veg for fi changes. if fi has increased,
! incorporate newly buried intercepted h2o into bottom-layer 
! snow, giving associated heat increment to soil, and mix the
! specific heat of newly buried veg (at tl) into tlsub. if fi
! has decreased, change temp of newly exhumed veg to tl, giving
! assoc heat increment to soil, and smear out intercepted h2o
      if (npn.gt.0) then
         do 600 i=kpti,kptj
            dfi(i) = fi(i) - fiold(i)
            if (dfi(i).gt.0.) then

! factor of xl*chl has been divided out of next line
               tlsub(i) = (tlsub(i)*fiold(i) + tl(i)*dfi(i)) / fi(i)
               zheat(i) = zheat(i) + dfi(i)*xl(i) * (                    &
                             wliql(i) * (                                &
                                ch2o*(tl(i)-tmelt)+hfus+cice*(           &
                                   tmelt-tsno(i,nsnolay)                 &
                                )                                        &
                             ) + wsnol(i) * cice*(tl(i)-tsno(i,nsnolay)) &
                          )

               hsno(i,nsnolay) = hsno(i,nsnolay)+dfi(i)*xl(i)*(wliql(i)+wsnol(i))/(rhos*fi(i))
            endif
            if (dfi(i).lt.0.) then
               zheat(i) = zheat(i) - dfi(i)*xl(i)*chl*(tlsub(i)-tl(i))
               rwork = (1.-fiold(i)) / (1.-fi(i))
               wliql(i) = wliql(i) * rwork
               wsnol(i) = wsnol(i) * rwork
            endif
600      continue
      endif

! areally average fluxes to be used by soil model. (don't use
! index due to mix call, but only need at all if npn > 0)
      if (npn.gt.0) then
         rwork = 1. / dtime
         do 700 i=kpti,kptj
            rwork2 = 1. - fiold(i)
            heatg(i) = rwork2*heatg(i)+fiold(i)*sflo(i,nsnolay+2)+zheat(i)*rwork

            solg(i)  = rwork2 * solg(i)
            fvapg(i) = rwork2 * fvapg(i)
            x1(i)    = rwork2 * raing(i)
            x2(i)    = zmelt(i)*rwork
            x3(i)    = tmelt
700      continue
         call mix(raing(kpti),traing(kpti), x1(kpti),traing(kpti), &
                  x2(kpti),x3(kpti), vzero,vzero,npt)
      endif
      return
end subroutine snow
