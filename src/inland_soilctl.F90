#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine soilctl(kpti, kptj)
! ---------------------------------------------------------------------
! steps soil/seaice model through one timestep
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control
      use inland_comatm
      use inland_comhyd
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d
      use inland_comcrop

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti    ! index of 1st point of little vector
                      ! in big lpt vector
      integer kptj    ! index of last point of little vector

! local variables
      integer npt     ! number of points in little vector
      integer i, k
      real*8 zfrez, & ! factor decreasing runoff fraction for tsoi < tmelt
             zrunf, & ! fraction of rain that doesn't stay in puddle (runoff fraction)
             zmax,  & ! maximum value of zrunf
             zmin,  & ! minimum value of zrunf
             grun1, & ! temporary storage for grunof
             rwork, & ! 
             wipre, & ! storing variable
             zdpud, & ! used to compute transfer from puddle to infiltration
             cx,    & ! average specific heat for soil, water and ice
             chav,  & ! average specific heat for water and ice
             zwsoi      

      real*8 owsoi(lbeg:lend,nsoilay), & ! old value of wsoi
             otsoi(lbeg:lend,nsoilay), & ! old value of tsoi
             c0pud(lbeg:lend,nsoilay), &  ! layer heat capacity due to puddles (=0 except for top)
             c1pud(lbeg:lend,nsoilay), &  ! updated av. specifilayer heat capacity due to  puddle
             fhtop(lbeg:lend),   &  ! heat flux through soil surface (for soilheat)
             fsqueez(lbeg:lend), &  ! excess amount of water (soilh2o) 
             dh(lbeg:lend),      &  ! correction if water at tsoi < tmelt or ice at temp > tmelt
             dw(lbeg:lend),      &  ! correction if water at tsoi < tmelt or ice at temp > tmelt
             zporos(lbeg:lend),  &
             dtsw                   ! Kaiyuan Li for Green-Ampt, delta volumatirc water content 

      integer layerno  ! Kaiyuan Li for Green-Ampt, layer number of wetting front           
!-----------------------------------------------------------------------
!     call setarrib (c0pud(kpti,1), npt*nsoilay, 0.0)
!     call setarrib (c1pud(kpti,1), npt*nsoilay, 0.0)
      npt = kptj - kpti + 1
      do k = 1, nsoilay
         do i = kpti, kptj
            c0pud(i,k) = 0.0
            c1pud(i,k) = 0.0
         end do
      end do


! for soil, set soil infiltration rate fwtop (for 
! soilh2o) and upper heat flux fhtop (for soilheat)
!
! also step puddle model wpud, wipud
!
! procedure is:
!
!   (0) immediately transfer any excess puddle liq to runoff
!
!   (1) apportion raing btwn puddle liquid(wpud) or runoff(grunof)
!
!   (2) apportion evap/condens (fvapg) btwn infil rate(fwtop), soil
!       ice (wisoi(i,1)), puddle liq (wpud), or puddle ice (wipud)
!
!   (3) transfer some puddle liquid to fwtop
!
!   (4) compute upper heat flx fhtop: includes fwtop*ch2o*tsoi(i,1)
!       to be consistent with whflo in soilheat, and accounts for
!       changing rain temp from traing to tsoi(i,1) and runoff temp
!       from tsoi to max(tsoi(i,1),tmelt)

      do 100 i = kpti, kptj

! (0) immediately transfer any excess puddle liq to runoff
!
! the following runoff formulation could give rise to very
! small amounts of negative runoff
     if (isimagro .eq. 0) then
         grunof(i) = min(wpud(i), max(dble(0.), wpud(i) + wipud(i) - wpudmax)) / dtime
         wpud(i) = wpud(i) - grunof(i) * dtime
     endif
! (1) apportion sfc-level rain between puddle liquid and runoff
!
! linear dependence of runoff fraction on wpud+wipud assumes
! uniform statistical distribution of sub-grid puddle 
! capacities between 0 and wpudmax. runoff fraction is 
! reduced linearly for tsoi < tmelt (by zfrez) in which case
! any rain will increase wpud which will be frozen to wipud
! below
         zfrez = max (dble(0.), min (dble(1.), (tsoi(i,1) - tmelt + .5) * 2.))
! always have some minimal amount of runoff (3%) even if 
! puddles are dry or soil is cold, but don't allow more than
! a specified maximum (30%), since the rain must also be allowed
! to infiltrate (and more surface runoff is generated later in
! step (5) anyway). zmin was "tuned" based on drainage-to-total
! runoff ratios for the Trout Lake region (for wpudmax=200mm).
! zmax is based on the assumption that 70% of precip goes to ET
! (on average over the Upper Midwest); under saturated conditions,
! the remainder (30%) is assumed to go directly to surface runoff.
! Both zmin and zmax can be considered tunable, but it would
! probably be better to just adjust wpudmax.
!
      if (isimagro .gt. 0) then

         zmin = 0.03
!        CJK 1-16-03      zmax = 1.00
         zmax = 0.30   ! zmax change according to TET on 1-16-03 
         zrunf = zmin + (zmax - zmin) * zfrez * max (dble(0.), min (dble(1.), (wpud(i) + wipud(i)) / wpudmax))
         wpud(i) = wpud(i) + (1. - zrunf) * raing(i) * dtime
         grunof(i) = zrunf * raing(i)

      else

!        always have some minimal amount of runoff (0.10) even if 
!        puddles are dry or soil is cold

         zrunf = zfrez * max (dble(0.), min (dble(1.), (wpud(i) + wipud(i)) / wpudmax))
         wpud(i) = wpud(i) + (1. - zrunf) * raing(i) * dtime
         grunof(i) = grunof(i) + zrunf * raing(i)
      endif

! (2) apportion evaporation or condensation between 4 h2o stores:
         rwork = fvapg(i) * dtime
         if (fvapg(i).ge.0.) then

! evaporation: split according to qglif
            fwtop(i) =          - qglif(i,1)*fvapg(i)
            wpud(i)  = wpud(i)  - qglif(i,3)*rwork
            wipud(i) = wipud(i) - qglif(i,4)*rwork
            wipre = wisoi(i,1)
            wisoi(i,1) = max(dble(0.),wipre-qglif(i,2)*rwork/(rhow*poros(i,1) * &
                             hsoi(1)))
            if (1.-wisoi(i,1).gt.epsilon) &
               wsoi(i,1) = wsoi(i,1)*(1.-wipre)/(1.-wisoi(i,1))
         else

! condensation: give all to puddles (to avoid wsoi, wisoi > 1)
            fwtop(i) = 0.
            wpud(i) = wpud(i)  - (qglif(i,1)+qglif(i,3))*rwork
            wipud(i)= wipud(i) - (qglif(i,2)+qglif(i,4))*rwork
         endif

! (3) transfer some puddle liquid to infiltration; can lead
!     to small amounts of negative wpud (in soilh2o) due to
!     round-off error

      if(isimagro .gt. 0) zdpud = rhow * dtime * max (dble(0.), 1.-wisoi(i,1))**2 * hydraul(i,1)

! isinfilt - Infiltration Function - 0: according to Darcy (1856); 1: according to Green-Ampt (1911) (default 0)
      if(isimagro .eq.0)then
       if (isinfilt.eq.0) then

          zdpud = rhow * dtime * max (dble(0.), 1.-wisoi(i,1))**2 * hydraul(i,1)
 
       else

! Added by Kaiyuan Li for incorporation of Green-ampt infiltration calculate potential infiltration (actual infiltration is fwpud)
! The Green-Ampt equation adopted bolow is from Julien et al. 1995
! Water resources buttetin vol. 31, No. 3: 523 - 536, 1995
 
         call delta_sw(fwpudtot(i), dtsw, layerno, i)   ! calculate delta soil water content

! calculate potential infiltration zdpud
          zdpud = 0.5 * ((max (dble(0.), 1.-wisoi(i,1))**2*hydraul(i, 1) * 1000.0 * dtime - 2.0 * fwpudtot(i)) &
                  +  sqrt((max (dble(0.), 1.-wisoi(i,1))**2*hydraul(i, 1) * 1000.0 * dtime - 2.0 * fwpudtot(i)) ** 2        &
                  +  8.0 * max (dble(0.), 1.-wisoi(i,1))**2*hydraul(i, 1) * 1000.0 * dtime * (cpwf(i, layerno) * 1000 *     &
                  dtsw + fwpudtot(i))))
       endif 
      endif !end of isimagro           


         fwpud(i) = max (dble(0.), min (wpud(i), zdpud)) / dtime
         c0pud(i,1) = ch2o*wpud(i) + cice*wipud(i)

! (4) compute upper soil heat flux
         fhtop(i) = heatg(i) + raing(i)*ch2o*(traing(i)-tsoi(i,1)) &
                  - grunof(i)*ch2o*max(tmelt-tsoi(i,1), dble(0.))
         soihfl(i) = fhtop(i)


! update diagnostic variables
         gadjust(i) = 0.0
100   continue

! reduce soil moisture due to transpiration (upsoi[u,l], from
! turvap).need to do that before other time stepping below since 
! specific heat of this transport is neglected
!
! first set porosflo, reduced porosity due to ice content, used
! as the effective porosity for uptake here and liquid hydraulics
! later in soilh2o. to avoid divide-by-zeros, use small epsilon
! limit; this will always cancel with epsilon or 0 in numerators
!
! also increment soil temperature to balance transpired water
! differential between temps of soil and leaf. physically
! should apply this to the tree, but would be awkward in turvap.
! 
! also, save old soil moisture owsoi and temperatures otsoi so
! implicit soilh2o and soilheat can aposteriori deduce fluxes.

     do 120 k = 1, nsoilay
         do 130 i = kpti, kptj
            porosflo(i,k) = poros(i,k) * max (epsilon, (1.-wisoi(i,k)))

! next line just for ice whose poros(i,k) is 0.0
            porosflo(i,k) = max (porosflo(i,k), epsilon)
            wsoi(i,k) = wsoi(i,k) - dtime * (upsoiu(i,k)+upsoil(i,k)) / &
                                            (rhow*porosflo(i,k)*hsoi(k))
            cx = c0pud(i,k) + hsoi(k) * ( &
                    (1.-poros(i,k))*csoi(i,k)*rhosoi(i,k)            &
                    + poros(i,k)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o*rhow &
                    + poros(i,k)*wisoi(i,k)*cice*rhow)

            tsoi(i,k) = tsoi(i,k) - dtime*ch2o*(upsoiu(i,k)*(tu(i)-tsoi(i,k)) &
                      + upsoil(i,k)*(tl(i)-tsoi(i,k))) / cx
            owsoi(i,k)  = wsoi(i,k)
            otsoi(i,k)  = tsoi(i,k)
130      continue
120   continue

! step soil moisture calculations
      call soilh2o (owsoi, fsqueez, kpti, kptj)
      
! set wsoi, wisoi to exactly 0 or 1 if differ by negligible 
! amount (needed to avoid epsilon errors in loop 400 below)
      call wadjust(kpti, kptj)

! update drainage and puddle
      do 200 i = kpti, kptj
        gdrain(i)  = wflo(i,nsoilay+1)
        c1pud(i,1) = ch2o*wpud(i) + cice*wipud(i)

! Kai
! isinfilt - Infiltration Function - 0: according to Darcy (1856); 1: according to Green-Ampt (1911) (default 0)
!
        if (isinfilt.eq.1) then    
       	   if ( raing(i) .lt. 0.00000001/3600.0)  then  !it was 0.0001/3600.0
       	      fwpudtot(i) = 0
           else
              fwpudtot(i) = fwpudtot(i) + (fwpud(i) - fsqueez(i)) * dtime
           end if
        end if

200   continue

! step temperatures due to conductive heat transport
      call soilheat (otsoi, owsoi, c0pud, fhtop, c1pud, kpti, kptj)

! heat-conserving adjustment for liquid/ice below/above melt
! point. uses exactly the same logic as for intercepted veg h2o
! in steph2o2. we employ the fiction here that soil liquid and
! soil ice both have density rhow, to avoid "pot-hole"
! difficulties of expansion on freezing. this is done by 
! dividing all eqns through by rhow(*hsoi).
!
! the factor (1-wsoi(old))/(1-wisoi(new)) in the wsoi increments
! results simply from conservation of h2o mass; recall wsoi is
! liquid content relative to ice-reduced pore space.
      do 400 k = 1, nsoilay
         do 410 i = kpti, kptj

! next line is just to avoid divide-by-zero for ice with
! poros = 0
         zporos(i) = max (poros(i,k), epsilon)
         rwork = c1pud(i,k)/rhow/hsoi(k) &
                 + (1.-zporos(i))*csoi(i,k)*rhosoi(i,k)/rhow
         chav = rwork &
                + zporos(i)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o &
                + zporos(i)*wisoi(i,k)*cice

! if liquid exists below melt point, freeze some to ice
!
! (note that if tsoi>tmelt or wsoi=0, nothing changes.)
! (also note if resulting wisoi=1, either dw=0 and prev
! wisoi=1, or prev wsoi=1, so use of epsilon is ok.)
         zwsoi = min (dble(1.), wsoi(i,k))
         dh(i) = chav * (tmelt-tsoi(i,k))
         dw(i) = min ( zporos(i)*(1.-wisoi(i,k))*zwsoi, max (dble(0.),dh(i)/hfus) )
         wisoi(i,k) = wisoi(i,k) +  dw(i)/zporos(i)
         wsoi(i,k)  = wsoi(i,k)  - (dw(i)/zporos(i))*(1.-zwsoi) &
                                    / max (epsilon,1.-wisoi(i,k))
         chav = rwork &
                 + zporos(i)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o &
                 + zporos(i)*wisoi(i,k)*cice
         tsoi(i,k) = tmelt - (dh(i)-hfus*dw(i)) / chav

! if ice exists above melt point, melt some to liquid
!
! note that if tsoi<tmelt or wisoi=0, nothing changes
!
! also note if resulting wisoi=1, dw=0 and prev wisoi=1,
! so use of epsilon is ok
         dh(i) = chav * (tsoi(i,k) - tmelt)
         dw(i) = min ( zporos(i)*wisoi(i,k), max (dble(0.), dh(i)/hfus) )
         wisoi(i,k) = wisoi(i,k) -  dw(i)/zporos(i)
         wsoi(i,k)  = wsoi(i,k)  + (dw(i)/zporos(i)) &
                      * (1.-wsoi(i,k)) / max(epsilon,1.-wisoi(i,k))
          chav = rwork &
                 + zporos(i)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o &
                 + zporos(i)*wisoi(i,k)*cice
          tsoi(i,k) = tmelt + (dh(i)-hfus*dw(i)) / chav

! reset porosflo (although not used after this)
         porosflo(i,k) = zporos(i) * max (epsilon, 1.-wisoi(i,k))
410      continue
400   continue
 
! set wsoi, wisoi to exactly 0 or 1 if differ by negligible 
! amount (roundoff error in loop 400 above can produce very
! small negative amounts)
      call wadjust(kpti, kptj)

! repeat ice/liquid adjustment for upper-layer puddles (don't 
! divide through by rhow*hsoi). upper-layer moistures wsoi,wisoi
! are already consistent with tsoi(i,1) > or < tmelt, and will 
! remain consistent here since tsoi(i,1) will not cross tmelt
      k = 1
      do 500 i = kpti, kptj

! if any puddle liquid below tmelt, freeze some to puddle ice
         rwork = ( (1.-poros(i,k))*csoi(i,k)*rhosoi(i,k) &
                 + poros(i,k)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o*rhow &
                 + poros(i,k)*wisoi(i,k)*cice*rhow) * hsoi(k)
         chav = ch2o*wpud(i) + cice*wipud(i) + rwork
         dh(i) = chav * (tmelt-tsoi(i,k))
         dw(i) = min (wpud(i), max (dble(0.), dh(i)/hfus))
         wipud(i) = wipud(i) + dw(i)
         wpud(i)  = wpud(i)  - dw(i)
         chav = ch2o*wpud(i) + cice*wipud(i) + rwork
         tsoi(i,k) = tmelt - (dh(i)-hfus*dw(i)) / chav
!
! (5) transfer any excess puddle liq to runoff
!
! the following runoff formulation could give rise to very
! small amounts of negative runoff
!
      if (isimagro .gt. 0) then
         grun1 = (min (wpud(i), max (0., wpud(i) + wipud(i) - wpudmax))) / dtime
!
         grunof(i) = grunof(i) + grun1
!
         wpud(i) = wpud(i) - grun1 * dtime
       endif
!
! if any puddle ice above tmelt, melt it and send to puddle liquid
! (not apportioned between puddle and surface runoff, to avoid
! potential double shunting to runoff, i.e. duplicating step 1).

! if any puddle ice above tmelt, melt some to puddle liquid
         dh(i) = chav * (tsoi(i,k)-tmelt)
         dw(i) = min (wipud(i), max (dble(0.), dh(i)/hfus))
         wipud(i) = wipud(i) - dw(i)
         wpud(i)  = wpud(i)  + dw(i)
         chav = ch2o*wpud(i) + cice*wipud(i) + rwork
         tsoi(i,k) = tmelt + (dh(i)-hfus*dw(i)) / chav
500   continue
    if(isimagro .eq. 0)then
!     if (ipointout .eq. 1 .and. myid .eq. 0) then
      if (ipointout .eq. 1) then
         i = 1001
         if (i .ge. kpti .and. i .le. kptj) then
            write(201,'(11(f9.3,1x),3(f12.3,1x))') &
            poros(i,1)*hsoi(1)*1e3*(wsoi(i,1)*(1 - wisoi(i,1)) + wisoi(i,1)), &
            poros(i,2)*hsoi(2)*1e3*(wsoi(i,2)*(1 - wisoi(i,2)) + wisoi(i,2)), &
            poros(i,3)*hsoi(3)*1e3*(wsoi(i,3)*(1 - wisoi(i,3)) + wisoi(i,3)), &
            poros(i,4)*hsoi(4)*1e3*(wsoi(i,4)*(1 - wisoi(i,4)) + wisoi(i,4)), &
            poros(i,5)*hsoi(5)*1e3*(wsoi(i,5)*(1 - wisoi(i,5)) + wisoi(i,5)), &
            poros(i,6)*hsoi(6)*1e3*(wsoi(i,6)*(1 - wisoi(i,6)) + wisoi(i,6)), &
            raina(i)*1e6,snowa(i)*1e6,grunof(i)*1e6,gdrain(i)*1e6,wpud(i),    &
            fi(i)*rhos*hsno(i,1)*1e3,fi(i)*rhos*hsno(i,2)*1e3,                &
            fi(i)*rhos*hsno(i,3)*1e3
!        call flush(101)
            write(202,'(6f8.2)') (tsoi(i,k),k=1,nsoilay)
!        call flush(202)
         end if

         i = 662
         if (i .ge. kpti .and. i .le. kptj) then
            write(1101,'(11(f9.3,1x),3(f12.3,1x))') &
            poros(i,1)*hsoi(1)*1e3*(wsoi(i,1)*(1 - wisoi(i,1)) + wisoi(i,1)), &
            poros(i,2)*hsoi(2)*1e3*(wsoi(i,2)*(1 - wisoi(i,2)) + wisoi(i,2)), &
            poros(i,3)*hsoi(3)*1e3*(wsoi(i,3)*(1 - wisoi(i,3)) + wisoi(i,3)), &
            poros(i,4)*hsoi(4)*1e3*(wsoi(i,4)*(1 - wisoi(i,4)) + wisoi(i,4)), &
            poros(i,5)*hsoi(5)*1e3*(wsoi(i,5)*(1 - wisoi(i,5)) + wisoi(i,5)), &
            poros(i,6)*hsoi(6)*1e3*(wsoi(i,6)*(1 - wisoi(i,6)) + wisoi(i,6)), &
            raina(i)*1e6,snowa(i)*1e6,grunof(i)*1e6,gdrain(i)*1e6,wpud(i),    &
            fi(i)*rhos*hsno(i,1)*1e3,fi(i)*rhos*hsno(i,2)*1e3,                &
            fi(i)*rhos*hsno(i,3)*1e3
!           call flush(1101)
            write(1102,'(6f8.2)') (tsoi(i,k),k=1,nsoilay)
!           call flush(1102)
         end if

         i = 1142
         if (i .ge. kpti .and. i .le. kptj) then
            write(2101,'(11(f9.3,1x),3(f12.3,1x))') &
            poros(i,1)*hsoi(1)*1e3*(wsoi(i,1)*(1 - wisoi(i,1)) + wisoi(i,1)), &
            poros(i,2)*hsoi(2)*1e3*(wsoi(i,2)*(1 - wisoi(i,2)) + wisoi(i,2)), &
            poros(i,3)*hsoi(3)*1e3*(wsoi(i,3)*(1 - wisoi(i,3)) + wisoi(i,3)), &
            poros(i,4)*hsoi(4)*1e3*(wsoi(i,4)*(1 - wisoi(i,4)) + wisoi(i,4)), &
            poros(i,5)*hsoi(5)*1e3*(wsoi(i,5)*(1 - wisoi(i,5)) + wisoi(i,5)), &
            poros(i,6)*hsoi(6)*1e3*(wsoi(i,6)*(1 - wisoi(i,6)) + wisoi(i,6)), &
            raina(i)*1e6,snowa(i)*1e6,grunof(i)*1e6,gdrain(i)*1e6,wpud(i),    &
            fi(i)*rhos*hsno(i,1)*1e3,fi(i)*rhos*hsno(i,2)*1e3,                &
            fi(i)*rhos*hsno(i,3)*1e3
!           call flush(2101)
            write(2102,'(6f8.2)') (tsoi(i,k),k=1,nsoilay)
!           call flush(2102)
         end if
      endif
     endif !end of isimagro
      return 
end subroutine soilctl
