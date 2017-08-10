#include "inland_config.h"
!---------------------------------------------------------------------
subroutine twoset(omega,betad,betai,avmu,gdir,coszen,iv,ib,loopi)
!---------------------------------------------------------------------
!
! sets two-stream parameters, given single-element transmittance
! and reflectance, leaf orientation weights, and cosine of the
! zenith angle, then adjusts for amounts of intercepted snow
!
! the two-stream parameters omega,betad,betai are weighted 
! combinations of the "exact" values for the 3 orientations:
! all vertical, all horizontal, or all random (ie, spherical)
!
! the vertical, horizontal weights are in oriev,orieh (comveg)
!
! the "exact" expressions are as derived in my notes(8/6/91,p.6).
! note that values for omega*betad and omega*betai are calculated
! and then divided by the new omega, since those products are 
! actually used in twostrib. also those depend *linearly* on the
! single-element transmittances and reflectances tauveg, rhoveg,
! which are themselves linear weights of leaf and stem values 
!
! for random orientation, omega*betad depends on coszen according
! to the function in array tablemu
!
! the procedure is approximate since omega*beta[d,i] and gdir
! should depend non-linearly on the complete leaf-angle
! distribution. then we should also treat leaf and stem angle
! distributions separately, and allow for the cylindrical
! shape of stems (norman and jarvis, app.b; the expressions 
! below are appropriate for flat leaves)
      use inland_parameters
      use inland_comveg
      use inland_com1d
      use inland_comcrop, only: isimagro

      implicit none
! ---------------------------------------------------------------------
! input-output variables
      integer loopi ! index of the little vector (1, numlv)
      integer ib, & ! waveband number (1= visible, 2= near-IR)
              iv    ! 1 for lower, 2 for upper story params (supplied)


! all quantities are returned unless otherwise noted:
      real*8 omega(lbeg:lend),  &  ! fraction of intercepted radiation that is scattered
             betad(lbeg:lend),  &  ! fraction of scattered *direct* radiation that is
!                                  !  scattered into upwards hemisphere
             betai(lbeg:lend),  &  ! fraction of scattered downward *diffuse* radiation
!                                  ! that is scattered into upwards hemisphere (or fraction
!                                  ! of scattered upward diffuse rad. into downwards hemis)
             avmu(lbeg:lend),   &  ! average diffuse optical depth
             gdir(lbeg:lend),   &  ! average projected leaf area into solar direction
             coszen(lbeg:lend)     ! cosine of solar zenith angle (supplied)

! local variables
      integer j, i, ntmu, itab, nsolz
      real*8 zrho, ztau, orand, ztab, rwork, y, o, x, betadsno, betaisno
      
      real*8 otmp(lbeg:lend)

      parameter (ntmu=100)
      real*8 tablemu(ntmu+1), omegasno(nband)
      save tablemu, omegasno, betadsno, betaisno

      data tablemu / 0.5000, 0.4967, 0.4933, 0.4900, 0.4867, 0.4833, 0.4800, &
                     0.4767, 0.4733, 0.4700, 0.4667, 0.4633, 0.4600, 0.4567, &
                     0.4533, 0.4500, 0.4467, 0.4433, 0.4400, 0.4367, 0.4333, &
                     0.4300, 0.4267, 0.4233, 0.4200, 0.4167, 0.4133, 0.4100, &
                     0.4067, 0.4033, 0.4000, 0.3967, 0.3933, 0.3900, 0.3867, &
                     0.3833, 0.3800, 0.3767, 0.3733, 0.3700, 0.3667, 0.3633, &
                     0.3600, 0.3567, 0.3533, 0.3500, 0.3467, 0.3433, 0.3400, &
                     0.3367, 0.3333, 0.3300, 0.3267, 0.3233, 0.3200, 0.3167, &
                     0.3133, 0.3100, 0.3067, 0.3033, 0.3000, 0.2967, 0.2933, &
                     0.2900, 0.2867, 0.2833, 0.2800, 0.2767, 0.2733, 0.2700, &
                     0.2667, 0.2633, 0.2600, 0.2567, 0.2533, 0.2500, 0.2467, &
                     0.2433, 0.2400, 0.2367, 0.2333, 0.2300, 0.2267, 0.2233, &
                     0.2200, 0.2167, 0.2133, 0.2100, 0.2067, 0.2033, 0.2000, &
                     0.1967, 0.1933, 0.1900, 0.1867, 0.1833, 0.1800, 0.1767, &
                     0.1733, 0.1700, 0.1667 /

      data omegasno /0.9, 0.7/
      data betadsno, betaisno /0.5, 0.5/
    
!    if(isimagro .gt. 0) then
!
! Assign leaf optical properties (taken from Sellers et al., 1996
! and Bonan, 1995)
! These are reflectance and transmission parameters depending on what part
! of the spectrum is used, what part of the canopy is used (lower or upper),
! and whether the leaves are green or brown
!
!      rhovegvlg = 0.10      ! vis leaf reflectance, lower story, green leaves
!      rhovegvlb = 0.36      ! vis leaf reflectance, lower story, brown leaves
!      rhovegvu = 0.10       ! vis leaf reflectance, upper story, green leaves
!
!      rhovegirlg = 0.48     ! nir leaf reflectance, lower story, green leaves
!      rhovegirlb = 0.58     ! nir leaf reflectance, lower story, brown leaves
!      rhovegiru = 0.40      ! nir leaf reflectance, upper story, green leaves
!
!      tauvegvlg = 0.07      ! vis leaf transmittance, lower story, green leaves
!      tauvegvlb = 0.22      ! vis leaf transmittance, lower story, brown leaves
!      tauvegvu = 0.05       ! vis leaf transmittance, upper story, green leaves
!
!      tauvegirlg = 0.25     ! nir leaf transmittance, lower story, green leaves
!      tauvegirlb = 0.38     ! nir leaf transmittance, lower story, brown leaves
!      tauvegiru = 0.20      ! nir leaf transmittance, upper story, green leaves
!    endif
!
! set two-stream parameters omega, betad, betai, gdir and avmu
! as weights of those for 100% vert,horiz,random orientations
!
      nsolz = nsol(loopi)
    if(isimagro .eq. 0) then
      zrho = rhoveg(ib,iv)
      ztau = tauveg(ib,iv)
    endif

      do 100 j=1,nsolz
         i = indsol(loopi,j)
    
    if(isimagro .gt. 0) then
!
! The following determines zrho (reflectance of an average leaf) and
! ztau (transmittance of an average leaf) for location i.
! rhoveg and tauveg are given above for both upper and lower
! canopies and for visible and near infrared wavebands. This new
! routine adjusts those initialized values for the lower canopy
! depending on how much of the canopy is green.
! zrho and ztau will be
! weighted by greenfracl (the fraction of lower canopy that is green)
! to allow values to go from full green values to full brown values.
! Note that zrho for near infrared is the same for both green and
! brown leaves but the calculation is given for consistency.
!
! iv is 1 for lower canopy and 2 for upper canopy
! ib is 1 for visible wavebands and 2 for near infrared wavebands
!
        if (iv.eq.2) then
          if (ib.eq.1) then
!
! visible values for the upper canopy
!
            zrho = rhovegvu
            ztau = tauvegvu
          else
!
! ir values for the upper canopy
!
            zrho = rhovegiru
            ztau = tauvegiru
          endif
        else
          if (ib.eq.1) then
!
! visible values for the lower canopy, weighted by how much of
! canopy is green
!
            zrho = greenfracl(i) * rhovegvlg + &
                   rhovegvlb * (1. - greenfracl(i))
!
            ztau = greenfracl(i) * tauvegvlg + &
                   tauvegvlb * (1. - greenfracl(i))
!
          else
!
! ir values for the lower canopy, weighted by how much of
! canopy is green
!
            zrho = greenfracl(i) * rhovegirlg + &
                   rhovegirlb * (1. - greenfracl(i))
!
            ztau = greenfracl(i) * tauvegirlg + &
                   tauvegirlb * (1. - greenfracl(i))
!
          endif
        endif
!
    endif ! check for crop existence
! weight for random orientation is 1 - those for vert and horiz
         orand = 1. - oriev(iv) - orieh(iv)
         omega(i) = zrho + ztau

! ztab is transmittance coeff - for random-orientation omega*betad,
! given by tablemu as a function of coszen
         itab = nint (coszen(i)*ntmu + 1)
         ztab = tablemu(itab)
         rwork = 1./omega(i)
         betad(i) = (oriev(iv) * 0.5*(zrho + ztau) + orieh(iv) * zrho + &
                    orand * ((1.-ztab)*zrho + ztab*ztau)) * rwork
         betai(i) = (oriev(iv) * 0.5*(zrho + ztau) + orieh(iv) * zrho + &
                    orand * ((2./3.)*zrho + (1./3.)*ztau)) * rwork
         gdir(i) = oriev(iv) * (2./pi)*sqrt(max(dble(0.),1.-coszen(i)*coszen(i))) + &
                   orieh(iv) * coszen(i) + orand * 0.5
         avmu(i) = 1.
100   continue

! adjust omega, betad and betai for amounts of intercepted snow
! (omegasno decreases to .6 of cold values within 1 deg of tmelt)
      if (iv.eq.1) then

! lower story
         do 210 j=1,nsolz
            i = indsol(loopi,j)
            y = fwetl(i)*(1.-rliql(i))
            o = omegasno(ib)*(.6 + .4*max(dble(0.),min(dble(1.),(tmelt-tl(i))/1.0)))
            otmp(i)  = omega(i)
            rwork = y * o
            omega(i) =  (1-y)*otmp(i)          + rwork
            betad(i) = ((1-y)*otmp(i)*betad(i) + rwork*betadsno) / omega(i)  
            betai(i) = ((1-y)*otmp(i)*betai(i) + rwork*betaisno) / omega(i)  
210      continue
      else

! upper story
         do 220 j=1,nsolz
            i = indsol(loopi,j)
            x = lai(i,iv) / max (lai(i,iv)+sai(i,iv), epsilon)
            y = x * fwetu(i)*(1.-rliqu(i)) + (1-x) *fwets(i)*(1.-rliqs(i))
            o = (x * min (dble(1.), max (dble(.6), (tmelt-tu(i))/0.1)) + (1-x) * &
                min (dble(1.), max (dble(.6), (tmelt-ts(i))/0.1))) * omegasno(ib) 

            otmp(i)  = omega(i)
            rwork = y * o
            omega(i) =  (1-y)*otmp(i)          + rwork
            betad(i) = ((1-y)*otmp(i)*betad(i) + rwork*betadsno) / omega(i)
            betai(i) = ((1-y)*otmp(i)*betai(i) + rwork*betaisno) / omega(i)
220      continue
      endif
      return
end subroutine twoset
