#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine setsoi(kpti, kptj)
! ---------------------------------------------------------------------
! sets diagnostic soil quantities
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comveg
      use inland_comsatp
      use inland_comcrop, only: isimagro

      implicit none
! ---------------------------------------------------------------------
! input variables
      integer kpti            ! index of 1st point of little vector 
!                             ! in big lpt vector
      integer kptj            ! index of last point of little vector

! local variables
      integer k, i, msand, mclay
      real*8 fsand,   &         ! fraction of sand in grid point
             fsilt,   &         ! fraction of silt in grid point
             fclay,   &         ! fraction of clay in grid point
             powliq,  &         ! liquid water content in fraction of soil depth
             powice,  &         ! ice water content in fraction of soil depth
             zcondry, &          ! dry-soil conductivity
             zwsoi,   &         ! volumetric water content of top soil layer 
             zvap,    &         ! latent heat of vaporisation at soil temp
             zsub,    &         ! latent heat of sublimation at soil temp
             zwpud,   &         ! fraction of soil surface covered by puddle
             rwork1, rwork2

! ---------------------------------------------------------------------
#define SETSOI_COMSAT
#include "inland_comsat.h"
! ---------------------------------------------------------------------
! set soil layer quantities
      do 100 k = 1, nsoilay
         do 110 i = kpti, kptj

! Convert input sand and clay percents to fractions
        if(isimagro .eq.0) then
! Convert input sand and clay percents to fractions
            msand = nint(sand(i,k))
            mclay = nint(clay(i,k))
            fsand = 0.01 * msand
            fclay = 0.01 * mclay
            fsilt = 0.01 * (100 - msand - mclay)
!
        endif
! update thermal conductivity (w m-1 k-1)
!
! based on c = c1**v1 * c2**v2 * c3**v3 * c4**v4 where c1,c2..
! are conductivities of soil grains, air, liquid and ice
! respectively, and v1,v2... are their volume fractions 
! (so v1 = 1-p where p is the porosity, and v1+v2+v3+v4 = 1).
! then condry = c1**(1-p) * c2**p  is the dry-soil
! conductivity, and c = condry * (c3/c2)**v3 * (c4/c2)**v4, 
! where c2 = conductivity of air = .025 w m-1 k-1.
! however this formula agrees better with williams+smith
! table 4 for wet (unfrozen) sand and clay if c2 is decreased
! to ~.005. (for peat in next section, ok if c2 = .025).
! also see lachenbruch etal,1982,jgr,87,9301 and refs therein.
            powliq = poros(i,k) * wsoi(i,k) * (1. - wisoi(i,k))
            powice = poros(i,k) * wisoi(i,k)

          if(isimagro .eq. 0) then
            zcondry = fsand * 0.300 + fsilt * 0.265 + fclay * 0.250
            consoi(i,k) = zcondry * ((0.56*100.)**powliq) * ((2.24*100.)**powice)
          else
              zcondry = fracsand(i,k) * 0.300 + fracsilt(i,k) * 0.265 + &
                        fracclay(i,k) * 0.250
!                       + forganic * 0.026      ! for future use CJK
              consoi(i,k) = zcondry * ((0.56*100.)**powliq) &
                        * ((2.24*100.)**powice)
          endif


!

110      continue
100   continue


! set qglif - the fraction of soil sfc evaporation from soil
! liquid (relative to total from liquid and ice)
!
! Changes by John Lenters.
! the partitioning of soil evaporation between puddle and soil 
! is now proportional to the (assumed) fraction of surface area covered by
! soil/puddle and liquid/ice. The parameter "zwpmax" is currently
! set to 0.5, but can be considered tunable. It is the assumed fraction of
! soil surface area covered by puddle when wpud+wipud>=wpudmax.

! zwpud:   fraction of surface area covered by puddle (range: 0 - zwpmax)
! zwpmax:  maximum value of zwpud (currently assumed to be 0.5)
! 1-zwpud: fraction of surface area covered by soil (range: (1-zwpmax) - 1.0)
! zwsoi:   volumetric water content of top soil layer (range: 0 - 1.0)
!
! qglif(i,1): fraction of soil evap (fvapg) from soil liquid
! qglif(i,2): fraction of soil evap (fvapg) from soil ice
! qglif(i,3): fraction of soil evap (fvapg) from puddle liquid
! qglif(i,4): fraction of soil evap (fvapg) from puddle ice
      do 200 i = kpti, kptj        
         zwpud = max (dble(0.0), min (zwpmax, zwpmax*(wpud(i)+wipud(i))/wpudmax) )
         zwsoi = (1. - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1)
         if (zwsoi.ge.epsilon) then

! for a wet surface
            rwork1 = 1.0 / zwsoi
            if (zwpud.ge.epsilon) then
               rwork2 = 1./(wpud(i) + wipud(i))
               qglif(i,1) = (1.-zwpud) * (1.-wisoi(i,1)) * wsoi(i,1) * rwork1
               qglif(i,2) = (1. - zwpud) * wisoi(i,1) * rwork1
               qglif(i,3) = zwpud * wpud(i) * rwork2
               qglif(i,4) = zwpud * wipud(i) * rwork2
            else
               qglif(i,1) = (1. - wisoi(i,1)) * wsoi(i,1) * rwork1
               qglif(i,2) = wisoi(i,1) * rwork1
               qglif(i,3) = 0.0
               qglif(i,4) = 0.0
            endif
         else

! for a 100% dry soil surface, assign all soil evap to the puddles.
! Note that for small puddle sizes, this could lead to negative
! puddle depths. However, for a 100% dry soil with small puddles,
! evaporation is likely to be very small or less than zero
! (condensation), so negative puddle depths are not likely to occur.
         if (zwpud.ge.epsilon) then
            rwork2 = 1./(wpud(i) + wipud(i))
            qglif(i,1) = 0.0
            qglif(i,2) = 0.0
            qglif(i,3) = zwpud * wpud(i) * rwork2
            qglif(i,4) = zwpud * wipud(i) * rwork2
         else
            if (tsoi(i,1).ge.tmelt) then

! above freezing
               qglif(i,1) = 0.
               qglif(i,2) = 0.
               qglif(i,3) = 1.
               qglif(i,4) = 0.
            else

! below freezing
               qglif(i,1) = 0.
               qglif(i,2) = 0.
               qglif(i,3) = 0.
               qglif(i,4) = 1.
            endif
         endif
      endif

! set latent heat values
      zvap = hvapf (tsoi(i,1), ta(i))
      zsub = hsubf (tsoi(i,1), ta(i))
      hvasug(i) = (qglif(i,1) + qglif(i,3)) * zvap + &
                  (qglif(i,2) + qglif(i,4)) * zsub 
      hvasui(i) = hsubf(tsno(i,1),ta(i))
200   continue   
      return
end subroutine setsoi
