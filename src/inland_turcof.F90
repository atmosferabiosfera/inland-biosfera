#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine turcof (iter, kpti, kptj)
! ---------------------------------------------------------------------
!
! solves for wind speeds at various levels
!
! also computes upper and lower-region air-air transfer coefficients
! and saves them in com1d arrays cu and cl for use by turvap,
! and similarly for the solid-air transfer coefficients
! su, ss, sl, sg and si
!
! iter = current iteration number
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d
      use inland_comcrop

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer iter             ! current iteration number
      integer kpti             ! index of 1st point of little vector
                               ! in big lpt vector
      integer kptj             ! index of last point of little vector
      integer jday

! local variables
      integer i

      real*8 x,tauu,a,b,c,d,taul,cai, cbi, cci,cdi, cei, cfi, &
             rwork,  &         ! working variable
             cdmax,  &         ! max value for cd
             ca,     &         ! to compute inverse air-air transfer coeffs
             sg0,    &         ! to compute air-solid transfer coeff for soil
             si0               ! to compute air-solid transfer coeff for ice     
      real*8 yu(kpti:kptj), yl(kpti:kptj), xfaca(lbeg:lend)

! set stratification factors for lower and upper regions
! using values from the previous iteration
      xfaca(kpti:kptj) = 1.0

      call fstrat(t34,t12,xfaca,q34,q12,z3,z2,alogl,alogl,alog2,u2,richl, &
                  straml,strahl,iter,kpti,kptj)
      call fstrat(t12,ta,tfaca,q12,qa,z1,za,alogu,alogu,aloga,ua,richu,stramu, &
                  strahu,iter,kpti,kptj)

! eliminate c/d from eq (28), tau_l/rho from (26),(27), to get
! lower-story roughness alogl. yl/bdl is (tau_l/rho)/(c+d)
!
! equation numbers correspond to lsx description section 4.e
      do 100 i = kpti, kptj
         x = ((alog4(i)-alogav(i))/vonk)**2 * bdl(i)
         rwork = 1. / expl(i)
         yl(i) = ((x+1)*expl(i) + (x-1)*rwork) / ((x+1)*expl(i) - (x-1)*rwork)
         alogl(i) = alog3(i) - vonk * sqrt(yl(i)/bdl(i))
100   continue 

! eliminate tau_l/rho from (24),(25), tau_u/rho and a/b from
! (22),(23), to get upper-story roughness alogu
! 
! yu/bdu is (tau_u/rho)/(a+b)
      do 110 i = kpti, kptj
         x = ((alog2(i)-alogl(i))/vonk)**2 * bdu(i) / straml(i)
         rwork = 1. / expu(i)
         yu(i) = ((x+1)*expu(i) + (x-1)*rwork) / ((x+1)*expu(i) - (x-1)*rwork)
         alogu(i) = alog1(i) - vonk * sqrt(yu(i)/bdu(i))
 110  continue

! define the maximum value of cd
      cdmax = 300.0 / (2.0 * dtime)

! get tauu (=tau_u/rho) from (21), a and b from (22),(23),
! taul (=tau_u/rho) from (25), c and d from (26),(27)
!
! changed the following to eliminate small errors associated with
! moving this code to single precision - affected c and d,
! which made u_ become undefined, as well as affecting some
! other variables
      do 200 i = kpti, kptj
         tauu = (ua(i) * vonk/(aloga(i)-alogu(i)))**2 * stramu(i)
         a = 0.5 * tauu * (yu(i)+1)/bdu(i)
         b = 0.5 * tauu * (yu(i)-1)/bdu(i)
         taul = bdu(i) * (a/expu(i) - b*expu(i))
         c = 0.5 * taul * (yl(i)+1)/bdl(i)
         d = 0.5 * taul * (yl(i)-1)/bdl(i)

! evaluate wind speeds at various levels, keeping a minimum 
! wind speed of 0.01 m/s at all levels


        if (isimagro .eq. 0) then
           u1(i)  = max (dble(0.01), sqrt (max (dble(0.0), (a+b))))
           u12(i) = max (dble(0.01), sqrt (max (dble(0.0), (a/exphu(i)+b*exphu(i)))))
           u2(i)  = max (dble(0.01), sqrt (max (dble(0.0), (a/expu(i) +b*expu(i)))))
           u3(i)  = max (dble(0.01), sqrt (max (dble(0.0), (c+d))))
           u34(i) = max (dble(0.01), sqrt (max (dble(0.0), (c/exphl(i)+d*exphl(i)))))
           u4(i)  = max (dble(0.01), sqrt (max (dble(0.0), (c/expl(i) +d*expl(i)))))
        else
           u1(i)  = max (0.01,ua(i)*0.98 )
           u12(i) = max (0.01,ua(i)*0.78 )
           u2(i)  = max (0.01,ua(i)*0.76 )
           u3(i)  = max (0.01,ua(i)*0.74 )
           u34(i) = max (0.01,ua(i)*0.72 )
           u4(i)  = max (0.01,ua(i)*0.64 )
        endif
 200  continue

!
! compute inverse air-air transfer coeffs
!
! use of inverse individual coeffs cai, cbi, cci, cdi, cei, cfi avoids
! divide-by-zero as vegetation vanishes - combine into
! upper-region coeff cu from za to z12, and lower-region coeff
! cl from z34 to z12, and also coeffs
      do 300 i = kpti, kptj
         ca = ua(i)*strahu(i)*vonk**2/((aloga(i)-alogu(i))*(aloga(i)-alog1(i)))
         ca = min (cdmax, ca / (1. + ca * 1.0e-20))
         cai = 1.0 / (rhoa(i)*ca)
         cbi = diu(i) * (z1(i)-z12(i)) / (rhoa(i) * 0.5*(u1(i)+u12(i)))
         cci = diu(i) * (z12(i)-z2(i)) / (rhoa(i) * 0.5*(u12(i)+u2(i)))
         cdi = (alog2(i)-alogl(i)) * (alog2(i)-alog3(i)) / &
               (rhoa(i)*u2(i)*strahl(i)*vonk**2)
         cei = dil(i) * (z3(i)-z34(i)) / (rhoa(i) * 0.5*(u3(i)+u34(i)))
         cfi = dil(i) * (z34(i)-z4(i)) / (rhoa(i) * 0.5*(u34(i)+u4(i)))
         cu(i) = 1.0 / (cai + cbi)
         cl(i) = 1.0 / (cci + cdi + cei)

! compute air-solid transfer coeffs for upper leaves, upper
! stems, lower story (su,ss,sl)
         su(i) = rhoa(i) * cleaf  * sqrt (u12(i) / dleaf(2))
         ss(i) = rhoa(i) * cstem  * sqrt (u12(i) / dstem(2))
         sl(i) = rhoa(i) * cgrass * sqrt (u34(i) / dleaf(1))

! compute air-solid transfer coeffs for soil and snow (sg,si)
!
! old technique
!
!       sg0 = rhoa(i) * u4(i) * (vonk/(alog4(i)-alogg(i)))**2
!       si0 = rhoa(i) * u4(i) * (vonk/(alog4(i)-alogi(i)))**2
!
! replace above formulations which depend on the log-wind profile
! (which may not work well below a canopy), with empirical formulation
! of Norman's. In the original LSX, turcof.f solves for the winds at
! the various levels from the momentum equations. This gives the transfer
! coefficients for heat and moisture. Heat and moisture eqns are then solved 
! in subroutine turvap. Using the empirical formulation of John Norman is 
! not consistent with the earlier solution for u4 (based on a logarithmic 
! profile just above the ground. However, this is used here because it 
! improved a lot simulations of the sensible heat flux over the 
! HAPEX-MOBILHY and FIFE sites
         sg0 = rhoa(i) * (0.004 + 0.012 * u4(i))
         si0 = rhoa(i) * (0.003 + 0.010 * u4(i))

! modify the cofficient to deal with cfi (see above)
         sg(i) = 1.0 / (cfi + 1.0 / sg0)
         si(i) = 1.0 / (cfi + 1.0 / si0)
300   continue
! JAF:  not necessary 
! if no veg, recalculate coefficients appropriately for a
! single logarithmic profile, and 2 fictitious levels just
! above soil/snow surface. these levels are arbitrary but are
! taken as z2 and z4, preset in vegdat to a few cm height
! for bare ground and ice. use strahu from above, which used
! t12 and alogu (ok after first iteration)
!
!     do 600 i = kpti, kptj
!
!       if ((fu(i).eq.0.0).and.(fl(i).eq.0.0)) then
!
!         z = rhoa(i)*ua(i)*strahu(i)*vonk**2 / (aloga(i)-alogav(i))
!
!         ca    = z / (aloga(i)-alog2(i))
!         cu(i) = rhoa(i)*min (cdmax,
!    >                          ca / (1. + ca / 1.0e+20))
!
!         cl(i) = z / (alog2(i)-alog4(i))
!
!         sg(i) = z / (alog4(i)-alogg(i))
!         si(i) = z / (alog4(i)-alogi(i))
!
!         alogu(i) = alogav(i)
!
!       endif
!
! 600 continue
      return
end subroutine turcof 
