#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine canini(kpti, kptj)
! ---------------------------------------------------------------------
!
! initializes aerodynamic quantities that remain constant 
! through one timestep
!
! note that some quantities actually are
! constant as long as the vegetation amounts and fractional
! coverage remain unchanged, so could re-arrange code for
! efficiency - currently all arrays initialized here are in
! com1d which can be overwritten elsewhere
!
! rwork is used throughout as a scratch variable to reduce number of
! computations
!
! ---------------------------------------------------------------------
! TODO: specify only the needed 'use' for every module
      use inland_com1d
#ifndef SINGLE_POINT_MODEL
      use inland_parameters, only: cair, cappa, cvap, grav, rair, rvap, dispuhf, &
                                ialoglhf, ialoguhf
#else /* SINGLE_POINT_MODEL */
      use inland_parameters, only: cair, cappa, cvap, rair, rvap, dispuhf,       &
                                ialoglhf, ialoguhf
      use inland_control, only: jday, istep
#endif /* SINGLE_POINT_MODEL */
      use inland_comveg
      use inland_comsoi
      use inland_comatm
      use inland_comsno

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector 
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

!-------------------------------Local variables-------------------------
      integer i

      real*8 bvegl, & ! e-folding depth in canopy for lower canopy
             cvegl, &
             dvegl, & ! diffusion coefficient for lower canopy
             bvegu, & ! e-folding depth in canopy for upper canopy 
             cvegu, &
             dvegu, & ! diffusion coefficient for upper canopy
             pa,    & ! pressure at level of atmospheric data
             siga,  & ! sigma level of atmospheric data
             rwork, & ! difference between top and bottom of canopy
             x,     & ! density of vegetation (without distinction betweenlai,sai) 
             x1         ! density of vegetation (different max)
!-----------------------------------------------------------------------
!
! define sigma level of atmospheric data
!   We are using once again SAGE's form to calculate: We assume a sigma value
! and based on this, we calculate the value of pa in the lowest atmosphere layer
! said, 10 meters above surface level.
!   For that height, since sage, the value is 0.999.
      siga = 0.999
!
! TODO: FIXME: transform tfaca(i) into constant tfac!
      tfaca(:) = 1.0 / (siga**cappa)

! atmospheric conditions at za
      do 100 i = kpti, kptj
         pa = psurf(i) * siga
         rhoa(i) = pa / (rair*ta(i)* (1.0 + (rvap/rair-1.0) * qa(i)))
         cp(i) = cair * (1.0 + (cvap/cair-1.0) * qa(i))
#ifndef SINGLE_POINT_MODEL
         za(i) = (psurf(i)-pa) / (rhoa(i)*grav)

! make sure that atmospheric level is higher than canopy top
         za(i) = max (za(i),ztop(i,2)+1.0)
#else /* SINGLE_POINT_MODEL */
! FIXME: 0D model uses fixed za to denote the observation tower's height.
!       We shall set this parameter in input files not hardcoded. - fzm

! Imbuzeiro: Fixed the tower height value (za) according the LBA-DMIP
! protocol. This parameter varies according site and user specifications
!         za(i) = 64.
#endif
100   continue 

! aerodynamic coefficients for the lower story
!
! cvegl (drag coeff for momentum) is proportional, and dvegl
! (diffusion coeff for momentum) inversely proportional,
! to x = density of vegetation (without distinction between
! lai,sai and fl*(1-fi)) - x is not allowed to be exactly
! zero to avoid divide-by-zeros, and for x>1 dvegl is 
! proportional to 1/x**2 so that roughness length tends to
! zero as x tends to infinity
!
! also the top, bottom and displacement heights z3(i),z4(i),
! displ(i) tend to particular values as the density tends to
! zero, to give same results as equations for no veg at all.
      do 200 i = kpti,kptj
         x = fl(i) * (1.0-fi(i)) * 2.0 * (lai(i,1)+sai(i,1)) / alaiml
         x  = min(x,dble(3.0))
         x1 = min(x,dble(1.0))

         rwork = max(ztop(i,1)-zbot(i,1),dble(0.01))
         cvegl = (0.4/rwork) * max(dble(1.e-5),x)
         dvegl = (0.1*rwork) / max(dble(1.e-5),x,x**2)

! e-folding depth in canopy
         bvegl = sqrt(dble(2.0)*cvegl/dvegl )

! [(tau/rho)/u**2] for inf canopy
         bdl(i) = 0.5*bvegl*dvegl

! 1 / diffusion coefficient
         dil(i) = 1./dvegl
         rwork = (1.0-x1) * (max(z0soi(i),z0sno)+0.01) 
         z3(i) = x1*ztop(i,1)+rwork
         z4(i) = x1*zbot(i,1)+rwork
         z34(i) = 0.5*(z3(i)+z4(i))
         exphl(i) = exp(0.5*bvegl*(z3(i)-z4(i)))
         expl(i)  = exphl(i)**2
         displ(i) = x1*0.7*z3(i)
200   continue 

! aerodynamic coefficients for the upper story
! same comments as for lower story
      do 300 i = kpti, kptj
         x = fu(i)*2.0*(lai(i,2)+sai(i,2)) / alaimu
         x  = min(x,dble(3.0))
         x1 = min(x,dble(1.0))

         rwork = max(ztop(i,2)-zbot(i,2),dble(.01))
         cvegu = (0.4/rwork) * max(dble(1.e-5),x)
         dvegu = (0.1 * rwork) / max(dble(1.e-5),x,x**2)
         rwork = 1./dvegu
         bvegu  = sqrt(2.0*cvegu*rwork)
         bdu(i) = 0.5*bvegu*dvegu
         diu(i) = rwork

         rwork = (1.0-x1) * (z3(i)+0.01)
         z1(i) = x1*ztop(i,2)+rwork
         z2(i) = x1*zbot(i,2)+rwork
         z12(i) = 0.5 * (z1(i) + z2(i))

         exphu(i) = exp(0.5 * bvegu * (z1(i) - z2(i)))
         expu(i)  = exphu(i)**2
         dispu(i) = x1 * dispuhf * z1(i) + (1.0 - x1) * displ(i)
300   continue 

! mixing-length logarithms
      do 400 i = kpti, kptj
         alogg(i)  = log(z0soi(i))
         alogi(i)  = log(z0sno)
         alogav(i) = (1.0 - fi(i)) * alogg(i) + fi(i) * alogi(i)

! alog4 must be > z0soi, z0sno to avoid possible problems later 
         alog4(i) = log(max(z4(i), 1.1*z0soi(i), 1.1*z0sno))

#ifdef SINGLE_POINT_MODEL
! alog3,2,1,a must have a logarythm
         if (jday.eq.1.and.istep.eq.1) then
            if (min(z3(i),z2(i)).le.displ(i)) then
               write(STDOUT,*) "*** WARNING: z3 or z2 (or both) are <= displ", &
                               " at ",__FILE__,":",__LINE__
                      
            endif
            if (min(z1(i),za(i)).le.dispu(i)) then
               write(STDOUT,*) "*** WARNING: z1 or za (or both) are <= dispu", &
                               " at ",__FILE__,":",__LINE__
            endif
         endif
#endif /* SINGLE_POINT_MODEL */

         alog3(i) = log(z3(i)-displ(i))
         alog2(i) = log(z2(i)-displ(i))
         alog1(i) = log(z1(i)-dispu(i))
         aloga(i) = log(za(i)-dispu(i))

! initialize u2, alogu, alogl for first iteration's fstrat
         u2(i)    = ua(i)/exphu(i)
         alogu(i) = log(max(dble(.01), ialoguhf*(z1(i)-z2(i))))
         alogl(i) = log(max(dble(.01), ialoglhf*(z3(i)-z4(i))))
400   continue
      return
end subroutine canini
