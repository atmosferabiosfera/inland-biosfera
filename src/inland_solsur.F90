#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine solsur(ib, loopi)
! ---------------------------------------------------------------------
! sets surface albedos for soil and snow, prior to other
! solar calculations
!
! ib = waveband number
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_com1d

      implicit none
! ---------------------------------------------------------------------
! input variable
      integer ib          ! waveband number. 1 = visible, 2 = near IR
      integer loopi       ! index of little vector in big vector

! local variable
      integer j, &    ! loop indice on number of points with >0 coszen
              nsolz, & ! number of points with >0 coszen
              i

      real*8 dinc,   & ! albedo correction du to soil moisture
             zw,     & ! liquid moisture content
             a7svlo, & ! snow albedo at low threshold temp., visible
             a7snlo, & ! snow albedo at low threshold temp., near IR
             a7svhi, & ! snow albedo at high threshold temp., visible
             a7snhi, & ! snow albedo at high threshold temp., near-IR
             t7shi,  & ! high threshold temperature for snow albedo
             t7slo     ! low  threshold temperature for snow albedo

      real*8 x(lbeg:lend), zfac(lbeg:lend)

! set the "standard" snow values:
      data a7svlo, a7svhi /0.90, 0.70/
      data a7snlo, a7snhi /0.60, 0.40/

!     t7shi ... high threshold temperature for snow albedo
!     t7slo ... low  threshold temperature for snow albedo
      t7shi = tmelt
      t7slo = tmelt - 15.0

! do nothing if all points in current strip have coszen le 0
      if (nsol(loopi).eq.0) then
         return
      endif

      nsolz = nsol(loopi)
      if (ib.eq.1) then
! soil albedos (visible waveband)
         do 100 j = 1, nsolz
            i = indsol(loopi,j)
! change the soil albedo as a function of soil moisture
            zw = wsoi(i,1) * (1.-wisoi(i,1))
            dinc = 1.0 + 1.0 * min (dble(1.), max (dble(0.0), 1. - (zw /.50) ))
            albsod(i) = min (albsav(i) * dinc, dble(.80))
            albsoi(i) = albsod(i)
100      continue

! snow albedos (visible waveband)
         do 110 j = 1, nsolz
            i = indsol(loopi,j)
            x(i) = (a7svhi*(tsno(i,1)-t7slo) + a7svlo*(t7shi-tsno(i,1))) / &
                   (t7shi-t7slo)
            x(i) = min (a7svlo, max (a7svhi, x(i)))
            zfac(i)   = max ( dble(0.), 1.5 / (1.0 + 4.*coszen(i)) - 0.5 )
            albsnd(i) = min ( dble(0.99), x(i) + (1.-x(i))*zfac(i))
            albsni(i) = min ( dble(1.), x(i))
110      continue
      else

! soil albedos (near-ir waveband)
         do 200 j = 1, nsolz
            i = indsol(loopi,j)

! LSX.2 formulation (different from lsx.1)
            zw = wsoi(i,1) * (1. - wisoi(i,1))
            dinc = 1.0 + 1.0 * min (dble(1.), max (dble(0.0), 1.0 - (zw / .50)  ))
            albsod(i) = min (albsan(i) * dinc, dble(.80))
            albsoi(i) = albsod(i)
200      continue

! snow albedos (near-ir waveband)
         do 210 j = 1, nsolz
            i = indsol(loopi,j)
            x(i) = (a7snhi * (tsno(i,1)-t7slo) + a7snlo * (t7shi - tsno(i,1))) &
                   / (t7shi - t7slo)
            x(i) = min (a7snlo, max (a7snhi, x(i)))
            zfac(i) = max ( dble(0.), 1.5/(1.+4.*coszen(i)) - 0.5 )
            albsnd(i) = min ( dble(0.99), x(i) + (1.-x(i))*zfac(i))
            albsni(i) = min ( dble(1.), x(i))
210      continue
      endif
      return
end subroutine solsur
