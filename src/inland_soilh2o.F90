#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine soilh2o(owsoi,fsqueez,kpti,kptj)
! ---------------------------------------------------------------------
!
! sets up call to tridiaib to solve implicit soil moisture eqn,
! using soil temperatures in wsoi (in comsoi)
!
! lower bc can be no h2o flow or free drainage, set by bperm below
!
! all arguments are supplied except wflo (returned):
!
! owsoi  = soil moistures at start of timestep
! fwtop  = h2o flux into top layer (infiltration)
! wflo   = downward h2o flow across layer boundaries
! itflo  = flow sub-timestep counter within overall surface step
! ntflo  = number of flow sub-timesteps per overall surface step
!
! local arrays and scalars:
!
! hsoim  = vertical distances between centers of layers
! wsoim  = interpolated moisture values at layer boundaries
! a,b,bwn[1] = intermediate terms (const for each pt -see notes)
! m,n    = exponents (constant for each point - see notes)
! e,f,g  = intermediate terms in algebraic devel (see notes)
! d1,2,3 = diagonals of tridiagonal systems of equations 
! rhs    = right-hand sides of systems of equations
! w1,2   = work arrays needed by tridiaib
! dmin   = minimum diffusivity for dry soils (m**2 s-1)
! rimp   = implicit fraction of the calculation (0 to 1)
! bperm  = 0 for impermeable (no drainage) lower bc,
!          1 for free drainage lower bc
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsoi
      use inland_com1d
      use inland_control
      use inland_comcrop, only:isimagro

      implicit none
!-----------------------------------------------------------------------
! input-output variables
!
! input variables
      integer kpti            ! index of 1st point of little vector
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

      real*8 owsoi(lbeg:lend,nsoilay), &! soil moistures at start of timestep
             fsqueez(lbeg:lend)         ! excess water at end of time step in soil column  

! local variables
      integer npt            ! number of points in little vector
      integer k, i, km1, kka, kkb , j
      real*8 dmin, &         ! minimum diffusivity for dry soils (m**2 s-1)
             rimp, &         ! implicit fraction of the calculation (0 to 1)
             zbex, z, dt, zz

      integer m(lbeg:lend), n(lbeg:lend) ! exponents 

      real*8 hsoim(nsoilay+1) ! vertical distances between centers of layers

      real*8 wsoim(lbeg:lend,nsoilay+1),&! interpolated moisture values at layer boundaries
             wsoia(lbeg:lend,nsoilay+1), wsoib(lbeg:lend,nsoilay+1), &
             weim(lbeg:lend,nsoilay+1), weip(lbeg:lend,nsoilay+1), &
             a(lbeg:lend),            & ! intermediate terms (const for each pt)
             b(lbeg:lend),bwn(lbeg:lend),bwn1(lbeg:lend), &
             e(lbeg:lend,nsoilay+1),  & ! intermediate terms in algebraic devel 
             f(lbeg:lend,nsoilay+1), g(lbeg:lend,nsoilay+1), &
             d1(mpt,nsoilay),         & ! diagonals of tridiagonal systems of equations 
             d2(mpt,nsoilay),d3(mpt,nsoilay), &
             rhs(mpt,nsoilay),        & ! right-hand sides of systems of equations
             wsoiw(mpt,nsoilay),      &
             w1(mpt,nsoilay),         & ! work arrays needed by tridia
             w2(mpt,nsoilay)
 
      save dmin, rimp
      data dmin, rimp /1.e-9, 1.0/
 
      npt = kptj - kpti + 1

! set level vertical distances, interpolated moistures, and
! interpolation weights
! top layer
      k = 1
      do 100 i = kpti, kptj
         hsoim(k) = 0.5 * hsoi(k)
         weim(i,k) = 0.0
         weip(i,k) = 1.0
         wsoim(i,k) = wsoi(i,k)
         wsoia(i,k) = min (wsoim(i,k), dble(1.0))
         wsoib(i,k) = min (wsoim(i,k), dble(1.0))
100   continue

! middle layers
      do 110 k = 2, nsoilay
         do 120 i = kpti, kptj
            hsoim(k) = 0.5 * (hsoi(k-1) + hsoi(k))
            weim(i,k) = 0.5 * hsoi(k) / hsoim(k)
            weip(i,k) = 1.0 - weim(i,k)
            wsoim(i,k) = weim(i,k) * wsoi(i,k-1) + weip(i,k) * wsoi(i,k)
            wsoia(i,k) = min (wsoim(i,k), dble(1.0))
            wsoib(i,k) = min (wsoim(i,k), dble(1.0))
120      continue
110   continue

! bottom layer
      k = nsoilay + 1
      do 130 i = kpti, kptj
         hsoim(k) = 0.5 * hsoi(k-1)
         weim(i,k) = 1.0
         weip(i,k) = 0.0
         wsoim(i,k) = wsoi(i,k-1)
         wsoia(i,k) = min (wsoim(i,k), dble(1.0))
         wsoib(i,k) = min (wsoim(i,k), dble(1.0))
130   continue

! set intermediate quantities e,f,g. these are terms in the
! expressions for the fluxes at boundaries between layers,
! so are zero for k=1. use bwn1 to account for minimum 
! diffusivity dmin. bperm is used for k=nsoilay+1 to set the
! type of the lower bc.
!
! top layer
      k = 1
      do i = kpti, kptj
         e(i,k) = 0.0
         f(i,k) = 0.0
         g(i,k) = 0.0
      end do

! middle layers
      do 200 k = 2, nsoilay
         do 210 i = kpti, kptj
! now that hydraul, suction and ibex can vary with depth,
! use averages of surrounding mid-layer values
!
! this is not rigorous since basic equations are in terms of
! potentials, not moisture - but that would require rewrite of soilh2o
!
! (see notes 8/27/93)
            a(i) = weim(i,k) * hydraul(i,k-1) + weip(i,k) * hydraul(i,k)
            b(i) = weim(i,k) * hydraul(i,k-1) * suction(i,k-1) * bex(i,k-1) + &
                   weip(i,k) * hydraul(i,k)   * suction(i,k)   * bex(i,k)
            zbex = weim(i,k) * bex(i,k-1) + weip(i,k) * bex(i,k) 
            m(i) = 2 * nint(zbex) + 3
            n(i) =     nint(zbex) + 2
            bwn1(i) = b(i) * (wsoib(i,k)**(n(i)-1))
            bwn(i)  = bwn1(i) * wsoib(i,k)

            if (bwn(i).lt.dmin) bwn1(i) = 0.0
            bwn(i) = max (bwn(i), dmin)
            e(i,k) = (-1.+rimp*m(i))*a(i)*(wsoia(i,k)**m(i)) + (        &
                        (1.-rimp)*bwn(i) - rimp*n(i)*bwn1(i)*wsoib(i,k) &
                     ) * (wsoi(i,k)-wsoi(i,k-1)) / hsoim(k)

            f(i,k) = - rimp*m(i)*a(i)*(wsoia(i,k)**(m(i)-1)) +          &
                     rimp*n(i)*bwn1(i)*(wsoi(i,k)-wsoi(i,k-1)) / hsoim(k)
            g(i,k) = rimp*bwn(i)
210      continue
200   continue

! bottom layer
      k = nsoilay + 1
      do 220 i = kpti, kptj
         a(i) = hydraul(i,nsoilay) 
         b(i) = hydraul(i,nsoilay)*suction(i,nsoilay)*ibex(i,nsoilay)
         m(i) = 2*ibex(i,nsoilay) + 3
         n(i) = ibex(i,nsoilay)   + 2
         e(i,k) = -a(i)*(wsoia(i,k)**m(i))*bperm
         f(i,k) = 0.0
         g(i,k) = 0.0
220   continue

! deduce all e,f,g in proportion to the minimum of the two 
! adjacent layers' (1-wisoi), to account for restriction of flow
! by soil ice. this will cancel in loop 300  with the factor 
! 1-wisoi in (one of) the layer's porosflo, even if wisoi=1 by 
! the use of epsilon limit. so a layer with wisoi=1 will form a 
! barrier to flow of liquid, but still have a predicted wsoi
      do 230 k = 1, nsoilay+1
         kka = max(k-1,1)
         kkb = min(k,nsoilay)
         do 240 i= kpti, kptj

! multiply by an additional factor of 1-wisoi for stability
            z = max(dble(0.),dble(1.)-max(wisoi(i,kka),wisoi(i,kkb)))**2

            e(i,k) = z * e(i,k)
            f(i,k) = z * f(i,k)
            g(i,k) = z * g(i,k)
240      continue
230   continue

! set matrix diagonals and right-hand sides
      do 300 k = 1, nsoilay
         do 310 i = kpti, kptj
            j = i - kpti + 1

            dt = dtime / (porosflo(i,k)*hsoi(k))
            d1(j,k) = dt*(f(i,k)*0.5*hsoi(k)/hsoim(k) - g(i,k)/hsoim(k))
            rhs(j,k) = wsoi(i,k) + dt*(e(i,k+1) - e(i,k))
310      continue
         if (k.eq.1) then
            do 320 i= kpti, kptj
               j = i - kpti + 1
               dt = dtime / (porosflo(i,k)*hsoi(k))
               rhs(j,k) = rhs(j,k) + dt*(fwtop(i)+fwpud(i))/rhow
320         continue
         endif

         if (k.lt.nsoilay) then
            km1 = max (k-1,1)
            do 330 i= kpti, kptj
               j = i - kpti + 1
               dt = dtime / (porosflo(i,k)*hsoi(k))
               d2(j,k) = 1. + dt*( &
                            - f(i,k+1)*0.5*hsoi(k+1)/hsoim(k+1) + f(i,k) *0.5* &
                            hsoi(km1)/hsoim(k) + g(i,k+1)/hsoim(k+1) + &
                            g(i,k) /hsoim(k) )
               d3(j,k) = dt*(-f(i,k+1)*0.5*hsoi(k)/hsoim(k+1)-g(i,k+1)/hsoim(k+1))
330         continue
         elseif (k.eq.nsoilay) then
            do 340 i= kpti, kptj
               j = i - kpti + 1
               dt = dtime / (porosflo(i,k)*hsoi(k))
               d2(j,k) = 1. + dt*(-f(i,k+1)+f(i,k)*0.5*hsoi(k-1)/hsoim(k)+g(i,k)/hsoim(k))
               d3(j,k) = 0.0
340         continue
         endif
300   continue

! solve the systems of equations

      call tridiaib(npt, npt, nsoilay, d1, d2, d3, rhs, wsoiw, w1, w2)

! Translates wsoiw(1 --> npt) from tridiaib to kpti --> kptj 
      do k = 1, nsoilay
         do i = kpti, kptj
            j = i - kpti + 1
            wsoi(i,k) = wsoiw(j,k)
         end do
      end do

      do 400 i = kpti, kptj
         fsqueez(i) = 0.0
         wflo(i,nsoilay+1) = - rhow * e(i,nsoilay+1)
400   continue

      do 500 k = nsoilay, 1, -1
         do 510 i = kpti, kptj
            zz = rhow * poros(i,k) * max(epsilon, (1.-wisoi(i,k))) * hsoi(k)
            wsoi(i,k) = wsoi(i,k) + dtime * fsqueez(i) / zz 
            fsqueez(i) = max (wsoi(i,k)-dble(1.),dble(0.)) * zz / dtime
            wsoi(i,k) = min (wsoi(i,k),dble(1.))
            wflo(i,k) = wflo(i,k+1) + (wsoi(i,k)-owsoi(i,k)) * zz / dtime
510      continue
500   continue

! step puddle liquid due to fsqueez and fwpud
!
! also subtract net puddle-to-top-layer flux from wflo(i,1),
! since puddle and top soil layer are lumped together in soilheat
! so upper wflo should be external flux only (evap/condens)
      do 600 i= kpti, kptj
         wpud(i)   = wpud(i)   + (fsqueez(i) - fwpud(i)) * dtime
         wflo(i,1) = wflo(i,1) + (fsqueez(i) - fwpud(i))
600   continue
      return
end subroutine soilh2o
