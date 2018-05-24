#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine irrad(kpti, kptj)
! ---------------------------------------------------------------------
! calculates overall emitted ir flux, and net absorbed minus
! emitted ir fluxes for upper leaves, upper stems, lower story,
! soil and snow. assumes upper leaves, upper stems and lower
! story each form a semi-transparent plane, with the upper-leaf
! plane just above the upper-stem plane. the soil and snow 
! surfaces have emissivities of 0.95
!
! the incoming flux is supplied in comatm array fira
!
! the emitted ir flux by overall surface system is returned in
! com1d array firb - the ir fluxes absorbed by upper leaves,
! upper stems, lower veg, soil and snow are returned in com1d 
! arrays firu, firs, firl, firg and firi
! 
! other com1d arrays used are:
!
! emu, ems, eml  = emissivities of the vegetation planes
! fup, fdown     = upward and downward fluxes below tree level
! ---------------------------------------------------------------------
      use inland_parameters, only: stef, avmuir
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d

      implicit none
!-----------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector 
                              ! in big lpt vector
      integer kptj            ! index of last point of little vector

! set emissivities of soil and snow
      real*8 emisoil, emisnow
      data emisoil, emisnow &
          /0.95,    0.95/

!     local arrays:
!     emu   = ir emissivity of upper-leaves veg plane
!     ems   = ir emissivity of upper-stems veg plane
!     eml   = ir emissivity of lower-story veg plane
!     emg   = ir emissivity (gray) of soil surface
!     emi   = ir emissivity (gray) of snow surface
!     fdown = downward ir flux below tree level per overall area
!     fup   = upward   ir flux below tree level per overall area
!     fdowng= downward ir flux below lower-story veg
!     fupg  = upward   ir flux below lower-story veg
!     fupgb = upward   ir flux above bare soil surface
!     fupi  = upward   ir flux above snow surface
      real*8 emu(kpti:kptj), ems(kpti:kptj), eml(kpti:kptj), emg(kpti:kptj), &
             emi(kpti:kptj), fdown(kpti:kptj), fdowng(kpti:kptj),            &
             fup(kpti:kptj), fupg(kpti:kptj), fupgb(kpti:kptj), fupi(kpti:kptj)

! use uniform value 1.0 for average diffuse optical depth
! (although an array for solar, all values are set to 1 in twoset).
      integer i
      do 100 i = kpti, kptj
         emu(i) = 1. - exp ( -lai(i,2) / avmuir )
         ems(i) = 1. - exp ( -sai(i,2) / avmuir )
         eml(i) = 1. - exp ( -(lai(i,1)+sai(i,1)) / avmuir )
         emg(i) = emisoil
         emi(i) = emisnow
         fdown(i) = (1.-fu(i)) * fira(i) + fu(i) * (        &
                       (1.-emu(i))*(1.-ems(i))*fira(i) +    &
                       emu(i)*(1.-ems(i))*stef*(tu(i)**4) + &
                       ems(i)*stef*(ts(i)**4)               &
                    )
         fdowng(i) = (1.-eml(i))*fdown(i)  + eml(i)*stef*(tl(i)**4)

         fupg(i)   = (1.-emg(i))*fdowng(i) + emg(i)*stef*(tg(i)**4)
         fupgb(i)  = (1.-emg(i))*fdown(i)  + emg(i)*stef*(tg(i)**4)
         fupi(i)   = (1.-emi(i))*fdown(i)  + emi(i)*stef*(ti(i)**4)
         
         fup(i) = (1.-fi(i))*(fl(i) *                                  &
                     (eml(i) *stef*(tl(i)**4) + (1.-eml(i))*fupg(i)) + &
                     (1.-fl(i))*fupgb(i)                               &
                  ) + fi(i) * fupi(i)

         firb(i) = (1.-fu(i)) * fup(i) + fu(i) * (                             &
                      (1.-emu(i))*(1.-ems(i))*fup(i) + emu(i)*stef*(tu(i)**4)+ &
                      ems(i)*(1.-emu(i))*stef*(ts(i)**4)                       &
                  )

         firu(i) = emu(i)*ems(i)*stef*(ts(i)**4) + emu(i)*(1.-ems(i))*fup(i) + &
                   emu(i)*fira(i) - 2*emu(i)*stef*(tu(i)**4)

         firs(i) = ems(i)*emu(i)*stef*(tu(i)**4) + ems(i)*fup(i) +       &
                   ems(i)*(1.-emu(i))*fira(i) - 2*ems(i)*stef*(ts(i)**4)

         firl(i) = eml(i)*fdown(i) + eml(i)*fupg(i) - 2*eml(i)*stef*(tl(i)**4)

         firg(i) = fl(i) * (fdowng(i) - fupg(i)) + (1.-fl(i)) * (fdown(i) - fupgb(i))

         firi(i) = fdown(i) - fupi(i)
100   continue
      return
end subroutine irrad
