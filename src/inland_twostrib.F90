#include "inland_config.h"
! ------------------------------------------------------------------------
subroutine twostrib(abvegd, abvegi, refld, refli, fbeldd, fbeldi, fbelid, &
                     fbelii, asurd, asuri, iv, coszen, ib, loopi)
! ------------------------------------------------------------------------
!
! solves canonical radiative transfer problem of two-stream veg
! layer + underlying surface of known albedo, for unit incoming
! direct or diffuse flux. returns flux absorbed within layer,
! reflected flux, and downward fluxes below layer. note that all
! direct fluxes are per unit horizontal zrea, ie, already 
! including a factor cos (zenith angle)
!
! the solutions for the twostream approximation follow Sellers (1985),
! and Bonan (1996) (the latter being the LSM documentation)
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comveg
      use inland_com1d

      implicit none
! ---------------------------------------------------------------------
!
! input-output variables
! [d,i] => per unit incoming direct, diffuse flux
!
      integer loopi           ! index of little vector in big vector
      integer ib
!
      real*8 abvegd(lbeg:lend),    & ! flux absorbed by two-stream layer (returned)
             abvegi(lbeg:lend),    & ! flux absorbed by two-stream layer (returned)
             refld(lbeg:lend),     & ! flux reflected above two-stream layer (returned)
             refli(lbeg:lend),     & ! flux reflected above two-stream layer (returned)
             fbeldd(lbeg:lend),    & ! downward direct  flux below two-stream layer(returned)
             fbeldi(lbeg:lend),    & ! downward direct  flux below two-stream layer(returned)
             fbelid(lbeg:lend),    & ! downward diffuse flux below two-stream layer(returned)
             fbelii(lbeg:lend),    & ! downward diffuse flux below two-stream layer(returned)
             asurd(lbeg:lend,nband),&! direct  albedo of underlying surface (supplied)
             asuri(lbeg:lend,nband),&! diffuse albedo of underlying surface (supplied)
             coszen(lbeg:lend)       ! cosine of direct zenith angle (supplied, must be gt 0)
 
!
! local variables
      integer i, j, iv, nsolz
      real*8 b, c, c0, d, f, h, k, q, p, sigma
      real*8 ud1, ui1, ud2, ui2, ud3, xai, s1, s2, p1, p2, p3, p4, rwork, dd1, &
             di1, dd2, di2, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, absurd,   &
             absuri
!
      real*8 omega(lbeg:lend),     & !
             betad(lbeg:lend),     & !
             betai(lbeg:lend),     & !
             avmu(lbeg:lend),      & !
             gdir(lbeg:lend),      & !
             tmp0(lbeg:lend)         !
!
! do nothing if all points in current strip have coszen le 0
      if (nsol(loopi).eq.0) return

      nsolz = nsol(loopi)
!
! calculate two-stream parameters omega, betad, betai, avmu, gdir
      call twoset (omega, betad, betai, avmu, gdir, coszen, iv, ib, loopi)

      do 100 j=1,nsolz
         i = indsol(loopi,j)
!
! the notations used here are taken from page 21 of Bonan's LSM documentation:
! Bonan, 1996: A Land Surface Model (LSM version 1.0) for ecological, hydrological,
! and atmospheric studies: Technical description and user's guide. NCAR Technical
! Note. NCAR/TN-417+STR, January 1996.
!
! some temporary variables are also introduced, which are from the original
! lsx model.
         b = 1. - omega(i) * (1.-betai(i))
         c = omega(i) * betai(i)
         tmp0(i) = b*b-c*c
         q = sqrt (max(dble(0.0), tmp0(i)) )
         k = gdir(i) / max(coszen(i), dble(0.01))
         p = avmu(i) * k
!
! next line perturbs p if p = q
         if ( abs(p-q) .lt. .001*p ) p = (1.+sign(dble(.001),dble(p-q))) * p

         c0 = omega(i) * p
         d = c0 * betad(i)
         f = c0 * (1.-betad(i))
         h = q / avmu(i)

         sigma = p*p - tmp0(i)
!
! direct & diffuse parameters are separately calculated
         ud1 = b - c/asurd(i,ib)
         ui1 = b - c/asuri(i,ib)
         ud2 = b - c*asurd(i,ib)
         ui2 = b - c*asuri(i,ib)
         ud3 = f + c*asurd(i,ib)

         xai = max (lai(i,iv) + sai(i,iv), epsilon)

         s1 = exp(-1.*h*xai)
         s2 = exp(-1.*k*xai)

         p1 = b + q
         p2 = b - q
         p3 = b + p
         p4 = b - p
         rwork = 1./s1
!
! direct & diffuse parameters are separately calculated
         dd1 = p1*(ud1-q)*rwork - p2*(ud1+q)*s1
         di1 = p1*(ui1-q)*rwork - p2*(ui1+q)*s1
         dd2 = (ud2+q)*rwork - (ud2-q)*s1
         di2 = (ui2+q)*rwork - (ui2-q)*s1
         h1 = -1.*d*p4 - c*f
         rwork = s2*(d-c-h1*(ud1+p)/sigma)
         h2 = 1./dd1*( (d-h1*p3/sigma)*(ud1-q)/s1 - p2*rwork )
         h3 = -1./dd1*( (d-h1*p3/sigma)*(ud1+q)*s1 - p1*rwork )
         h4 = -1.*f*p3 - c*d
         rwork = s2*(ud3-h4*(ud2-p)/sigma)
         h5 = -1./dd2*( h4*(ud2+q)/(sigma*s1) + rwork )
         h6 = 1./dd2*( h4*s1*(ud2-q)/sigma + rwork )
         h7 = c*(ui1-q)/(di1*s1)
         h8 = -1.*c*s1*(ui1+q)/di1
         h9 = (ui2+q)/(di2*s1)
         h10= -1.*s1*(ui2-q)/di2
!
! save downward direct, diffuse fluxes below two-stream layer
         fbeldd(i) = s2
         fbeldi(i) = 0.
         fbelid(i) = h4/sigma*s2 + h5*s1 + h6/s1
         fbelii(i) = h9*s1 + h10/s1
!
! save reflected flux, and flux absorbed by two-stream layer
!
         refld(i) = h1/sigma + h2 + h3
         refli(i) = h7 + h8
         absurd = (1.-asurd(i,ib)) * fbeldd(i) + (1.-asuri(i,ib)) * fbelid(i)
         absuri = (1.-asuri(i,ib)) * fbelii(i)
         abvegd(i) = max (dble(0.), dble(1.) - refld(i) - absurd)
         abvegi(i) = max (dble(0.), dble(1.) - refli(i) - absuri)
!
! if no veg, make sure abveg (flux absorbed by veg) is exactly zero
! if this is not done, roundoff error causes small (+/-)
! sols, soll values in solarf and subsequent problems in turvap
! via stomataib
         if (xai.lt.epsilon) abvegd(i) = 0.0
         if (xai.lt.epsilon) abvegi(i) = 0.0
!
! some terms needed in canopy scaling
! the canopy scaling algorithm assumes that the net photosynthesis
! is proportional to absored par (apar) during the daytime. during night,
! the respiration is scaled using a 10-day running-average daytime canopy
! scaling parameter.
!
! apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
!
! in the equations below, 
!
!   k = term[u,l](i,6)
!   h = term[u,l](i,7)
!
!   A = term[u,l](i,1) * ipardir(0)
!   B = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
!   C = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
!
! calculations performed only for visible (ib=1)
!
         if (ib.eq.1) then
            if (iv.eq.1) then
               terml(i,1) = k * (1. + (h4-h1) / sigma)
               terml(i,2) = h * (h5 - h2)
               terml(i,3) = h * (h9 - h7)
               terml(i,4) = h * (h3 - h6)
               terml(i,5) = h * (h8 - h10)
               terml(i,6) = k
               terml(i,7) = h
            else
               termu(i,1) = k * (1. + (h4-h1) / sigma)
               termu(i,2) = h * (h5 - h2)
               termu(i,3) = h * (h9 - h7)
               termu(i,4) = h * (h3 - h6)
               termu(i,5) = h * (h8 - h10)
               termu(i,6) = k
               termu(i,7) = h
            endif
         end if
100   continue
      return
end subroutine twostrib
