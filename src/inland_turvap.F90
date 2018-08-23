#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine turvap(iter,niter,kpti,kptj,xu,xs,xl,chux,chsx,chlx,chgx,wlgx,wigx, &
                  cog,coi,zirg,ziri,wu,ws,wl,wg,wi,tuold,tsold,tlold,tgold, &
                  tiold,qgfac0)
! ---------------------------------------------------------------------
!
! solves canopy system with linearized implicit sensible heat and
! moisture fluxes
!
! first, assembles matrix arr of coeffs in linearized equations
! for tu,ts,tl,t12,t34,q12,q34,tg,ti and assembles the right hand
! sides in the rhs vector
!
! then calls linsolve to solve this system, passing template mplate of
! zeros of arr
!
! finally calculates the implied fluxes and stores them
! for the agcm, soil, snow models and budget calcs
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control
      use inland_comatm
      use inland_comhyd
      use inland_comsatp
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d
      use inland_comcrop

      implicit none
!-----------------------------------------------------------------------
!
! input variables
      integer niter, iter
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer i,j, k
      integer npt      ! number of points in little vector

      real*8 rwork, zwtot, rwork2, tgav, tiav, tuav,&
             tsav, tlav, quav, qsav, qlav, qgav, qiav, zwpud, zwsoi,&
             psig, hfac, hfac2, zwopt, zwdry, betaw, emisoil, ee1, qs1,&
             dqs1, xnumer, xdenom, betafac, betas

      real*8 fradu(kpti:kptj),frads(kpti:kptj), fradl(kpti:kptj),qu(kpti:kptj),&
             qs(kpti:kptj),ql(kpti:kptj),qg(kpti:kptj),qi(kpti:kptj), &
             dqu(kpti:kptj),dqs(kpti:kptj),dql(kpti:kptj),dqg(kpti:kptj), &
             dqi(kpti:kptj),tupre(kpti:kptj),tspre(kpti:kptj),tlpre(kpti:kptj),&
             tgpre(kpti:kptj),tipre(kpti:kptj),suw(kpti:kptj),ssw(kpti:kptj), &
             slw(kpti:kptj),sut(kpti:kptj),slt(kpti:kptj),slt0(kpti:kptj), &
             suh(kpti:kptj),ssh(kpti:kptj),slh(kpti:kptj),qgfac(kpti:kptj)

! FIXME: these arrays should have dimension kpti:kptj instead of lbeg:lend ,
!        and then iterate over j = i - kpti + 1
      real*8 xu(lbeg:lend),xs(lbeg:lend),xl(lbeg:lend),chux(lbeg:lend), &
             chsx(lbeg:lend),chlx(lbeg:lend),chgx(lbeg:lend),wlgx(lbeg:lend), &
             wigx(lbeg:lend),cog(lbeg:lend),coi(lbeg:lend),zirg(lbeg:lend), &
             ziri(lbeg:lend),wu(lbeg:lend),ws(lbeg:lend),wl(lbeg:lend), &
             wg(lbeg:lend),wi(lbeg:lend),tuold(lbeg:lend),tsold(lbeg:lend), &
             tlold(lbeg:lend),tgold(lbeg:lend),tiold(lbeg:lend), &
             qgfac0(lbeg:lend)


      integer nqn
      parameter (nqn=9)
!      real*8 arr(mpt,nqn,nqn), &
!             rhs(mpt,nqn), &
!             vec(mpt,nqn)
      real*8 arr(kptj-kpti+1,nqn,nqn), &
             rhs(kptj-kpti+1,nqn), &
             vec(kptj-kpti+1,nqn)
      integer mplate(nqn,nqn)

!                  tu  ts  tl t12 t34 q12 q34  tg  ti
!                  ----------------------------------
      data mplate / 1,  0,  0,  1,  0,  1,  0,  0,  0, & !tu
                    0,  1,  0,  1,  0,  1,  0,  0,  0, & !ts
                    0,  0,  1,  0,  1,  0,  1,  0,  0, & !tl
                    1,  1,  0,  1,  1,  0,  0,  0,  0, & !t12
                    0,  0,  1,  1,  1,  0,  0,  1,  1, & !t34
                    1,  1,  0,  0,  0,  1,  1,  0,  0, & !q12
                    0,  0,  1,  0,  0,  1,  1,  1,  1, & !q34
                    0,  0,  0,  0,  1,  0,  1,  1,  0, & !tg
                    0,  0,  0,  0,  1,  0,  1,  0,  1  / !ti

!-----------------------------------------------------------------------
#include "inland_comsat.h"

!-----------------------------------------------------------------------
!
! if first iteration, save original canopy temps in t*old
! (can use tsoi,tsno for original soil/snow skin temps), for
! rhs heat capacity terms in matrix soln, and for adjustment
! of canopy temps after each iteration
!
! also initialize soil/snow skin temps tg, ti to top-layer temps
!
! the variables t12, t34, q12, q34, for the first iteration
! are saved via global arrays from the previous gcm timestep,
! this is worth doing only if the agcm forcing is
! smoothly varying from timestep to timestep
      npt = kptj - kpti + 1
      if (iter.eq.1) then

! weights for canopy coverages
      do 10 i = kpti, kptj
         xu(i) = 2.0 * lai(i,2) * fu(i)
         xs(i) = 2.0 * sai(i,2) * fu(i)
         xl(i) = 2.0 * (lai(i,1) + sai(i,1)) * fl(i) * (1.0 - fi(i))
10    continue

! specific heats per leaf/stem area
      do 20 i = kpti, kptj
         chux(i) = chu + ch2o * wliqu(i) + cice * wsnou(i)
         chsx(i) = chs + ch2o * wliqs(i) + cice * wsnos(i)
         chlx(i) = chl + ch2o * wliql(i) + cice * wsnol(i)
20    continue

      do 30 i = kpti, kptj
         rwork = poros(i,1) * rhow
         chgx(i) = ch2o * wpud(i) + cice * wipud(i) &
                   + ((1.-poros(i,1)) * csoi(i,1) * rhosoi(i,1) &
                   + rwork * (1.-wisoi(i,1)) * wsoi(i,1) * ch2o &
                   + rwork * wisoi(i,1) * cice) * hsoi(1)

         wlgx(i) = wpud(i) + rwork * (1. - wisoi(i,1)) * wsoi(i,1) * hsoi(1)

         wigx(i) = wipud(i) + rwork * wisoi(i,1) * hsoi(1)
30      continue

! conductivity coeffs between ground skin and first layer
         do 40 i = kpti, kptj
            cog(i) = consoi(i,1) / (0.5 * hsoi( 1))
            coi(i) = consno      / (0.5 * max (hsno(i,1), hsnotop))
40       continue

! d(ir emitted) / dt for soil
         rwork = 4. * 0.95 * stef
         do 50 i = kpti, kptj
            zirg(i) = rwork * (tg(i)**3)
            ziri(i) = rwork * (ti(i)**3)
50       continue

! updated temperature memory
         do 60 i = kpti, kptj
            tuold(i) = tu(i)
            tsold(i) = ts(i)
            tlold(i) = tl(i)
            tgold(i) = tg(i)
            tiold(i) = ti(i)
60       continue
      endif

! set implicit/explicit factors w* (0 to 1) for this iteration
! w* is 1 for fully implicit, 0 for fully explicit
! for first iteration, impexp and impexp2 set w* to 1
      call impexp (wu, tu, chux, wliqu, wsnou, iter, kpti, kptj)
      call impexp (ws, ts, chsx, wliqs, wsnos, iter, kpti, kptj)
      call impexp (wl, tl, chlx, wliql, wsnol, iter, kpti, kptj)
      call impexp (wg, tg, chgx, wlgx,  wigx,  iter, kpti, kptj)

! call impexp2 for snow model
      call impexp2 (wi, ti, tiold, iter, kpti, kptj)

! adjust t* for this iteration
!
! in this routine we are free to choose them,
! since they are just the central values about which the
! equations are linearized - heat is conserved in the matrix
! solution because t*old are used for the rhs heat capacities
!
! here, let t* represent the previous soln if it was fully
! implicit, but weight towards t*old depending on the amount
! (1-w*) the previous soln was explicit
!
! this weighting is necessary for melting/freezing surfaces, for which t*
! is kept at t*old, presumably at or near tmelt


      do 80 i = kpti, kptj
         tu(i) = wu(i) * tu(i) + (1.0 - wu(i)) * tuold(i)
         ts(i) = ws(i) * ts(i) + (1.0 - ws(i)) * tsold(i)
         tl(i) = wl(i) * tl(i) + (1.0 - wl(i)) * tlold(i)
         tg(i) = wg(i) * tg(i) + (1.0 - wg(i)) * tgold(i)
         ti(i) = wi(i) * ti(i) + (1.0 - wi(i)) * tiold(i)
80    continue

! save current "central" values for final flux calculations
      do 90 i = kpti, kptj
         tupre(i) = tu(i)
         tspre(i) = ts(i)
         tlpre(i) = tl(i)
         tgpre(i) = tg(i)
         tipre(i) = ti(i)
90    continue

! calculate various terms occurring in the linearized eqns,
! using values of t12, t34, q12, q34 from
! the previous iteration
!
! specific humidities for canopy and ground, and derivs wrt t
! for canopy
!
! limit derivs to avoid -ve implicit q's below,
! as long as d(temp)s in one iteration are le 10 deg k
      do 100 i = kpti, kptj
         ee1      = esat(tu(i))
         qu(i)  = qsat (ee1, psurf(i))
         dqu(i) = dqsat (tu(i), qu(i))
         dqu(i) = min (dqu(i), qu(i) * 0.1)

         ee1      = esat(ts(i))
         qs(i)  = qsat (ee1, psurf(i))
         dqs(i) = dqsat (ts(i), qs(i))
         dqs(i) = min (dqs(i), qs(i) * 0.1)

         ee1      = esat(tl(i))
         ql(i)  = qsat (ee1, psurf(i))
         dql(i) = dqsat (tl(i), ql(i))
         dql(i) = min (dql(i), ql(i) * 0.1)

         ee1      = esat(tg(i))
         qg(i)  = qsat (ee1, psurf(i))
         dqg(i) = dqsat (tg(i), qg(i))
         dqg(i) = min (dqg(i), qg(i) * 0.1)

         ee1      = esat(ti(i))
         qi(i)  = qsat (ee1, psurf(i))
         dqi(i) = dqsat (ti(i), qi(i))
         dqi(i) = min (dqi(i), qi(i) * 0.1)
100   continue

! set qgfac0, factor by which soil surface specific humidity
! is less than saturation
!
! it is important to note that the qgfac expression should
! satisfy timestep cfl criterion for upper-layer soil moisture
! for small wsoi(i,1)
!
! for each iteration, qgfac is set to qgfac0, or to 1 if
! condensation onto soil is anticipated (loop 110 in canopy.f)
!
! Evaporation from bare soil is calculated using the "beta method"
! (e.g., eqns 5 & 7 of Mahfouf and Noilhan 1991, JAM 30 1354-1365),
! but converted to the "alpha method" (eqns 2 & 3 of M&N), to match
! the structure in inland. The conversion from the beta to alpha
! method is through the relationship:
!   alpha * qgs - q34 = beta * (hfac * qgs - q34),
! from which one solves for alpha (which is equal to qgfac0):
!   qgfac0 = alpha = (beta * hfac) + (1 - beta)*(q34/qgs)
      do 105 i = kpti, kptj

! first calculate the total saturated fraction at the soil surface
! (including puddles ... see soil.f)
         zwpud = max (dble(0.0), min (dble(0.5), 0.5 * (wpud(i) + wipud(i))/wpudmax))
         zwsoi = (1.0 - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1)
         zwtot = zwpud + (1. - zwpud) * zwsoi

! next calculate the matric potential (from eqn 9.3 of Campbell and
! Norman), multiply by gravitational acceleration to get in units
! of J/kg, and calculate the relative humidity at the soil water
! surface (i.e., within the soil matrix), based on thermodynamic
! theory (eqn 4.13 of C&N)
         psig = -grav * suction(i,1) * (zwtot ** (-bex(i,1)))
         hfac = exp(psig/(rvap*tg(i)))

! then calculate the relative humidity of the air (relative to
! saturation at the soil temperature). Note that if hfac2 > 1
! (which would imply condensation), then qgfac is set to 1
! later in the code (to allow condensation to proceed at the
! "potential rate")
!gabriel apagar
if (i.eq.85) write(*,*) "hfac2",hfac2,"q34(i)",q34(i),"qg(i)",qg(i),1
         hfac2 = q34(i)/qg(i)
if (i.eq.85) write(*,*) "hfac2",hfac2,"q34(i)",q34(i),"qg(i)",qg(i),2

! set the "beta" factor and then calculate "alpha" (i.e., qgfac0)
! as the beta-weighted average of the soil water RH and the "air RH"
! First calculate beta_w:
         zwopt = 1.0
         zwdry = swilt(i,1)
         betaw = max(dble(0.0), min(dble(1.), (zwtot - zwdry)/(zwopt - zwdry)))
       if(isimagro .gt. 0)then
!
! limit evap if soil is frozen or snow-covered
          if (tg(i).le.273.16.or.fi(i).gt.0.0) then
            betaw = 0.01
          endif
       endif
! Next convert beta_w to beta_s (see Milly 1992, JClim 5 209-226):
         emisoil = 0.95
         ee1      = esat(t34(i))
         qs1    = qsat (ee1, psurf(i))
         dqs1   = dqsat (t34(i), qs1)
         xnumer = hvap * dqs1
         xdenom = cp(i) + (4.0 * emisoil * stef * (t34(i))**3) / sg(i)
         betafac = xnumer / xdenom
         betas = betaw / (1.0 + betafac * (1.0 - betaw))

! Combine hfac and hfac2 into qgfac0 ("alpha") using beta_s
!gabriel apagar
if (i.eq.85) write(*,*) "qgfac0(i)",qgfac0(i),"betas",betas,"hfac",hfac,"hfac2",hfac2,1
         qgfac0(i) = betas * hfac + (1. - betas) * hfac2
if (i.eq.85) write(*,*) "qgfac0(i)",qgfac0(i),"betas",betas,"hfac",hfac,"hfac2",hfac2,2
105   continue

! set fractions covered by intercepted h2o to 1 if dew forms
!
! these fwet*x are used only in turvap, and are distinct from
! the real fractions fwet* that are set in fwetcal
!
! they must be exactly 1 if q12 > qu or q34 > ql, to zero transpiration
! by the factor 1-fwet[u,l]x below, so preventing "-ve" transp
!
! similarly, set qgfac, allowing for anticipated dew formation
! to avoid excessive dew formation (which then infiltrates) onto
! dry soils
      do 110 i = kpti, kptj
         fwetux(i) = fwetu(i)
         if (q12(i).gt.qu(i)) fwetux(i) = 1.0
         fwetsx(i) = fwets(i)
         if (q12(i).gt.qs(i)) fwetsx(i) = 1.0
         fwetlx(i) = fwetl(i)
         if (q34(i).gt.ql(i)) fwetlx(i) = 1.0
         qgfac(i) = qgfac0(i)
         if (q34(i).gt.qg(i)) qgfac(i) = 1.0

! set net absorbed radiative fluxes for canopy components
         fradu(i) = 0.0
         if (lai(i,2).gt.epsilon) &
            fradu(i) = (solu(i) + firu(i)) / (2.0 * lai(i,2))
         frads(i) = 0.0
         if (sai(i,2).gt.epsilon) &
            frads(i) = (sols(i) + firs(i)) / (2.0 * sai(i,2))
         fradl(i) = 0.0
         if ((lai(i,1)+sai(i,1)).gt.epsilon) &
            fradl(i) = (soll(i) + firl(i)) / &
            (2.0 * (lai(i,1) + sai(i,1)))
110   continue

! calculate canopy-air moisture transfer coeffs for wetted
! leaf/stem areas, and for dry (transpiring) leaf areas
!
! the wetted-area coeffs suw,ssw,slw are constrained to be less
! than what would evaporate 0.8 * the intercepted h2o mass in
! this timestep (using previous iteration's q* values)
!
! this should virtually eliminate evaporation-overshoots and the need
! for the "negative intercepted h2o"  correction in steph2o2
      do 200 i = kpti, kptj

! coefficient for evaporation from wet surfaces in the upper canopy:
         suw(i) = min ( fwetux(i) * su(i), &
                  0.8 * (wliqu(i) + wsnou(i)) / &
                  max (dtime * (qu(i) - q12(i)), epsilon))

! coefficient for transpiration from average upper canopy leaves:
         sut(i) = (1.0 - fwetux(i)) * 0.5 * &
                  (totcondub(i) * frac(i,1) + &
                   totcondub(i) * frac(i,2) + &
                   totcondub(i) * frac(i,3) + &
                   totconduc(i) * frac(i,4) + &
                   totcondub(i) * frac(i,5) + &
                   totconduc(i) * frac(i,6) + &
                   totcondub(i) * frac(i,7) + &
                   totcondub(i) * frac(i,8))
         sut(i) = max (dble(0.0), sut(i))

! coefficient for sensible heat flux from upper canopy:
         suh(i) = suw(i) * (rliqu(i)  * hvapf(tu(i),ta(i)) + &
                  (1.-rliqu(i)) * hsubf(tu(i),ta(i))) + &
                  sut(i) * hvapf(tu(i),ta(i))

! coefficient for evaporation from wet surfaces on the stems:
         ssw(i) = min (fwetsx(i) * ss(i), &
                  0.8 * (wliqs(i) + wsnos(i)) &
                  / max (dtime * (qs(i) - q12(i)), epsilon))

! coefficient for sensible heat flux from stems:
         ssh(i) = ssw(i) * (rliqs(i)  * hvapf(ts(i),ta(i)) + &
                  (1.-rliqs(i)) * hsubf(ts(i),ta(i)))


! coefficient for evaporation from wet surfaces in the lower canopy:
!
! for the regional U.S. version, this is multiplied by the green
! fraction to reduce interception evaporation during the winter
! and spring, which otherwise seems too high. There is some
! physical basis for this, especially in snow covered regions
! where snow cover can cause masking and/or compaction of grass,
! thereby reducing interception. (This should be accounted for
! by fi(i) in ginvap (see below), but areal snow cover fraction
! is usually significantly underestimated in IBIS.)
!

       if (isimagro .gt. 0) then
           slw(i) = min (fwetlx(i) * sl(i) * min(0.1, greenfracl(i)), &
                    0.8 * (wliql(i) + wsnol(i)) &
                    / max (dtime * (ql(i) - q34(i)), epsilon))
! coefficient for transpiration from average lower canopy leaves:
           slt0(i) = (1. - fwetlx(i)) * 0.5 * &
                     (totcondls(i) * frac(i,9) + &
                     totcondls(i) * frac(i,10) + &
                     totcondl4(i) * frac(i,11) + &
                     totcondl3(i) * frac(i,12) + &
                     totcondc3(i) * frac(i,13) + &
                     totcondc4(i) * frac(i,14) + &
                     totcondc3(i) * frac(i,15) + &
                     totcondc4(i) * frac(i,16))
       else
           slw(i) = min (fwetlx(i) * sl(i), &
                    0.8 * (wliql(i) + wsnol(i)) &
                    / max (dtime * (ql(i) - q34(i)), epsilon))

!
! coefficient for transpiration from average lower canopy leaves:
           slt0(i) = (1. - fwetlx(i)) * 0.5 * &
                     (totcondls(i) * frac(i,9) + &
                     totcondls(i) * frac(i,10) + &
                     totcondl4(i) * frac(i,11) + &
                     totcondl3(i) * frac(i,12))

       endif ! check for crop existence

         slt0(i) = max (dble(0.), slt0(i))

! averaged over stems and lower canopy leaves:
         slt(i) = slt0(i) * lai(i,1) / max (lai(i,1)+sai(i,1), epsilon)

! coefficient for sensible heat flux from lower canopy:
         slh(i) = slw(i) * (  rliql(i)  * hvapf(tl(i),ta(i)) + &
                  (1.-rliql(i)) * hsubf(tl(i),ta(i))) + &
                  slt(i) * hvapf(tl(i),ta(i))
200   continue

! set the matrix of coefficients and the right-hand sides
! of the linearized equations.
! arr, rhs and vec are defined from 1 -> npt, instead of kpti --> kptj
! j is used to translate them.
      arr(:,:,:) = 0.0
      rhs(:,:) = 0.0
      rwork = 1. / dtime

! upper leaf temperature tu
      do 300 i = kpti, kptj
         j = i - kpti + 1
         rwork2 = su(i) * cp(i)
         arr(j,1,1) = chux(i) * rwork + wu(i) * rwork2 &
                      + wu(i) * suh(i) * dqu(i)
         arr(j,1,4) = -rwork2
         arr(j,1,6) = -suh(i)
         rhs(j,1) = tuold(i) * chux(i) * rwork &
                    - (1.-wu(i)) * rwork2 * tu(i) &
                    - suh(i) * (qu(i) - wu(i) * dqu(i) * tu(i)) &
                    + fradu(i) - pfluxu(i)
300   continue

! upper stem temperature ts
      do 310 i = kpti, kptj
         j = i - kpti + 1
         rwork2 = ss(i)*cp(i)
         arr(j,2,2) = chsx(i) * rwork &
                      + ws(i) * rwork2 &
                      + ws(i) * ssh(i) * dqs(i)
         arr(j,2,4) = -rwork2
         arr(j,2,6) = -ssh(i)
         rhs(j,2) = tsold(i) * chsx(i) * rwork &
                    - (1.-ws(i)) * rwork2 * ts(i) &
                    - ssh(i) * (qs(i)-ws(i) * dqs(i) * ts(i)) &
                    + frads(i) - pfluxs(i)
310   continue

! lower veg temperature tl
      do 320 i = kpti, kptj
         j = i - kpti + 1
         rwork2 = sl(i) * cp(i)
         arr(j,3,3) = chlx(i) * rwork + wl(i) * rwork2 + wl(i) * slh(i) * dql(i)
         arr(j,3,5) = -rwork2
         arr(j,3,7) = -slh(i)
         rhs(j,3) = tlold(i) * chlx(i) * rwork &
                    - (1.-wl(i)) * rwork2 * tl(i) &
                    - slh(i) * (ql(i) - wl(i) * dql(i) * tl(i)) &
                    + fradl(i) - pfluxl(i)
320   continue

! upper air temperature t12
      do 330 i = kpti, kptj
         j = i - kpti + 1
         rwork = xu(i) * su(i)
         rwork2 = xs(i) * ss(i)
         arr(j,4,1) = -wu(i) * rwork
         arr(j,4,2) = -ws(i) * rwork2
         arr(j,4,4) = cu(i) + cl(i) + rwork + rwork2
         arr(j,4,5) = -cl(i)
         rhs(j,4) = cu(i) * ta(i) * tfaca(i) + (1.-wu(i)) * rwork * tu(i) &
                    + (1.-ws(i)) * rwork2 * ts(i)
330   continue

! lower air temperature t34
      do 340 i = kpti, kptj
         j = i - kpti + 1
         rwork = xl(i)*sl(i)
         rwork2 = fi(i)*si(i)
         arr(j,5,3) = -wl(i)*rwork
         arr(j,5,4) = -cl(i)
         arr(j,5,5) = cl(i) + rwork + (1.-fi(i)) * sg(i) + rwork2
         arr(j,5,8) = -wg(i) * (1.-fi(i)) * sg(i)
         arr(j,5,9) = -wi(i) * rwork2
         rhs(j,5) = (1.-wl(i)) * rwork * tl(i) &
                    + (1.-wg(i)) * (1.-fi(i)) * sg(i) * tg(i) &
                    + (1.-wi(i)) * rwork2 * ti(i)
340   continue

! upper air specific humidity q12
      do 350 i = kpti, kptj
         j = i - kpti + 1
         rwork = xu(i) * (suw(i) + sut(i))
         rwork2 = xs(i) * ssw(i)
         arr(j,6,1) = -wu(i) * rwork * dqu(i)
         arr(j,6,2) = -ws(i) * rwork2 * dqs(i)
         arr(j,6,6) = cu(i) + cl(i) + rwork + rwork2
         arr(j,6,7) = -cl(i)
         rhs(j,6) = cu(i) * qa(i) &
                   + rwork  * (qu(i) - wu(i) * dqu(i) * tu(i)) &
                   + rwork2 * (qs(i) - ws(i) * dqs(i) * ts(i))
350   continue

! lower air specific humidity q34
      do 360 i = kpti, kptj
         j = i - kpti + 1
         rwork  = xl(i) * (slw(i) + slt(i))
         rwork2 = (1.- fi(i)) * sg(i)
         arr(j,7,3) = -wl(i) * rwork * dql(i)
         arr(j,7,6) = -cl(i)
         arr(j,7,7) = cl(i) + rwork + rwork2 + fi(i) * si(i)
         arr(j,7,8) = -wg(i) * rwork2 * qgfac(i) * dqg(i)
         arr(j,7,9) = -wi(i) * fi(i) * si(i) * dqi(i)
         rhs(j,7) = rwork * (ql(i) - wl(i) * dql(i) * tl(i)) &
                    + rwork2 * qgfac(i) * (qg(i) - wg(i) * dqg(i) * tg(i)) &
                    + fi(i) * si(i) * (qi(i) - wi(i) * dqi(i) * ti(i))
360   continue

! soil skin temperature
!
! (there is no wg in this eqn since it solves for a fully
! implicit tg. wg can be thought of as the fractional soil
! area using a fully implicit soln, and 1-wg as that using a
! fully explicit soln. the combined soil temperature is felt
! by the lower air, so wg occurs in the t34,q34 eqns above.)
      do 370 i = kpti, kptj
         j = i - kpti + 1
         rwork  = sg(i) * cp(i)
         rwork2 = sg(i) * hvasug(i)
         arr(j,8,5) = -rwork
         arr(j,8,7) = -rwork2
         arr(j,8,8) = rwork + rwork2 * qgfac(i) * dqg(i) &
                     + cog(i) + zirg(i)
         rhs(j,8) = -rwork2 * qgfac(i) * (qg(i) - dqg(i) * tg(i)) &
                   + cog(i) * tsoi(i,1) &
                   + solg(i) + firg(i) + zirg(i) * tgold(i)
370   continue

! snow skin temperature
!
! (there is no wi here, for the same reason as for wg above.)
      do 380 i = kpti, kptj
         j = i - kpti + 1
         rwork  = si(i) * cp(i)
         rwork2 = si(i) * hvasui(i)
         arr(j,9,5) = -rwork
         arr(j,9,7) = -rwork2
         arr(j,9,9) = rwork + rwork2 * dqi(i) + coi(i) + ziri(i)
         rhs(j,9) = -rwork2 * (qi(i) - dqi(i) * ti(i)) &
                    + coi(i) * tsno(i,1) &
                    + soli(i) + firi(i) + ziri(i) * tiold(i)
380   continue

! solve the systems of equations (on the little vector 1 --> npt)
      call linsolve (arr, rhs, vec, mplate, nqn, npt)

! translate this iteration's solution to t*, q12, q34 (1 --> npt to
! kpti --> kptj)
      do 400 i = kpti, kptj
         j = i - kpti + 1
         tu(i)  = vec(j,1)
         ts(i)  = vec(j,2)
         tl(i)  = vec(j,3)
         t12(i) = vec(j,4)
         t34(i) = vec(j,5)
         tg(i)  = vec(j,8)
         ti(i)  = vec(j,9)
         q12(i) = vec(j,6)
         q34(i) = vec(j,7)
!gabriel apagar
if (i.eq.85) then
         write(*,*) "tu(i) ",tu(i)
         write(*,*) "ts(i) ",ts(i)
         write(*,*) "tl(i) ",tl(i)
         write(*,*) "t12(i)",t12(i)
         write(*,*) "t34(i)",t34(i)
         write(*,*) "tg(i) ",tg(i)
         write(*,*) "ti(i) ",ti(i)
         write(*,*) "q12(i)",q12(i)
         write(*,*) "q34(i)",q34(i)
end if
400   con
nue

! all done except for final flux calculations,
! so loop back for the next iteration (except the last)
      if (iter.lt.niter) return

! evaluate sensible heat and moisture fluxes (per unit
! leaf/stem/snow-free/snow-covered area as appropriate)
!
! *******************************
! diagnostic sensible heat fluxes
! *******************************
      do 500 i = kpti, kptj
         fsena(i) = cp(i) * cu(i) * (ta(i) * tfaca(i) - t12(i))
         tgav = wg(i) * tg(i) + (1.- wg(i)) * tgpre(i)
         fseng(i) = cp(i) * sg(i) * (tgav - t34(i))
         tiav = wi(i) * ti(i) + (1.- wi(i)) * tipre(i)
         fseni(i) = cp(i) * si(i) * (tiav - t34(i))
         tuav = wu(i) * tu(i) + (1. - wu(i)) * tupre(i)
         fsenu(i) = cp(i) * su(i) * (tuav - t12(i))
         tsav = ws(i) * ts(i) + (1. - ws(i)) * tspre(i)
         fsens(i) = cp(i) * ss(i) * (tsav - t12(i))
         tlav = wl(i) * tl(i) + (1. - wl(i)) * tlpre(i)
         fsenl(i) = cp(i) * sl(i) * (tlav - t12(i))
500   continue

! *************************
! calculate moisture fluxes
! *************************
      do 510 i = kpti, kptj

! total evapotranspiration from the entire column
         fvapa(i)  = cu(i) * (qa(i)-q12(i))

! evaporation from wet surfaces in the upper canopy
! and transpiration per unit leaf area - upper canopy
         quav = qu(i) + wu(i) * dqu(i) * (tu(i) - tupre(i))
         fvapuw(i) = suw(i) * (quav-q12(i))
         fvaput(i) = max (dble(0.0), sut(i) * (quav-q12(i)))

! evaporation from wet surfaces on stems
         qsav = qs(i) + ws(i) * dqs(i) * (ts(i) - tspre(i))
         fvaps(i) = ssw(i) * (qsav - q12(i))

! evaporation from wet surfaces in the lower canopy
! and transpiration per unit leaf area - lower canopy
         qlav = ql(i) + wl(i) * dql(i) * (tl(i) - tlpre(i))
         fvaplw(i) = slw(i) * (qlav-q34(i))
         fvaplt(i) = max (dble(0.0), slt0(i) * (qlav-q34(i)))

! evaporation from the ground
         qgav = qg(i) + wg(i) * dqg(i) * (tg(i) - tgpre(i))
         !gabriel apagar
         if (i.eq.85) write(*,*) 'fvapg(i)',fvapg(i),'sg(i)',sg(i),'q34(i)',q34(i),'qgfac(i)',qgfac(i),'qgav',qgav,1
         fvapg(i) = sg(i) * (qgfac(i) * qgav - q34(i))
         if (i.eq.85) write(*,*) 'fvapg(i)',fvapg(i),'sg(i)',sg(i),'q34(i)',q34(i),'qgfac(i)',qgfac(i),'qgav',qgav,1


! evaporation from the snow
         qiav = qi(i) + wi(i) * dqi(i) * (ti(i) - tipre(i))
         fvapi(i) = si(i) * (qiav-q34(i))
510   continue

! adjust ir fluxes
      do 520 i = kpti, kptj
         firg(i) = firg(i) - wg(i)*zirg(i)*(tg(i) - tgold(i))
         firi(i) = firi(i) - wi(i)*ziri(i)*(ti(i) - tiold(i))
         firb(i) = firb(i) + (1.-fi(i))*wg(i)*zirg(i)*(tg(i) - tgold(i)) &
                   + fi(i) *wi(i)*ziri(i)*(ti(i) - tiold(i))

! impose constraint on skin temperature
         ti(i) = min (ti(i), tmelt)
520   continue

! set upsoi[u,l], the actual soil water uptake rates from each
! soil layer due to transpiration in the upper and lower stories,
! for the soil model
      do 600 k = 1, nsoilay
         do 610 i = kpti, kptj

! isimrwu - Root water uptake module - 0: according to Foley et al., 1996; 1 according to Li et al. (2006) (default 0)

         if (isimrwu .eq. 1) then

            upsoiu(i,k) = fvaput(i) * 2.0 * lai(i,2) * fu(i) *    &
                          stressu(i,k) / max (stre_tu(i), epsilon)

            upsoil(i,k) = fvaplt(i) * 2.0 * lai(i,1) * fl(i) *    &
                          (1. - fi(i)) *                           &
                          stressl(i,k) / max (stre_tl(i), epsilon)
         else

            upsoiu(i,k) = fvaput(i) * 2.0 * lai(i,2) * fu(i) * &
                          stressu(i,k) / max (stresstu(i), epsilon)
            upsoil(i,k) = fvaplt(i) * 2.0 * lai(i,1) * fl(i) * &
                          (1. - fi(i)) * &
                          stressl(i,k) / max (stresstl(i), epsilon)
         endif

610      continue
600   continue

! set net evaporation from intercepted water, net evaporation
! from the surface, and net transpiration rates
      do 700 i = kpti, kptj

! evaporation from intercepted water
         ginvap(i) = fvapuw(i) * 2.0 * lai(i,2) * fu(i) + &
                     fvaps (i) * 2.0 * sai(i,2) * fu(i) + &
                     fvaplw(i) * 2.0 * (lai(i,1) + sai(i,1)) * &
                     fl(i) * (1. - fi(i))

! evaporation from soil and snow surfaces
         gsuvap(i) = fvapg(i)  * (1. - fi(i)) + fvapi(i)  * fi(i)

! transpiration
         gtrans(i) = fvaput(i) * 2.0 * lai(i,2) * fu(i) + &
                     fvaplt(i) * 2.0 * lai(i,1) * fl(i) * (1.-fi(i))

         gtransu(i) = fvaput(i) * 2.0 * lai(i,2) * fu(i)
         gtransl(i) = fvaplt(i) * 2.0 * lai(i,1) * fl(i) * (1.-fi(i))
700   continue

!     if (ipointout .eq. 1 .and. myid .eq. 0) then
      if (ipointout .eq. 1) then
         i = 1001
         if (i .ge. kpti .and. i .le. kptj) then
            write(104,'(4f7.2)') gtransu(i)*1e6,gtransl(i)*1e6,gsuvap(i)* &
                                 1e6,ginvap(i)*1e6
!           call flush(104)
         end if
         i = 662
         if (i .ge. kpti .and. i .le. kptj) then
            write(1104,'(4f7.2)') gtransu(i)*1e6,gtransl(i)*1e6,gsuvap(i)* &
                                  1e6,ginvap(i)*1e6
!           call flush(1104)
         end if
         i = 1142
         if (i .ge. kpti .and. i .le. kptj) then
            write(2104,'(4f7.2)') gtransu(i)*1e6,gtransl(i)*1e6,gsuvap(i)* &
                                  1e6,ginvap(i)*1e6
!           call flush(2104)
         end if
      end if
      return
end subroutine turvap
