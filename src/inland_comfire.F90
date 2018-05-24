#include "inland_config.h"
module inland_comfire
      implicit none
      public
      save

!
! -------
! comfire
! -------
!

      real*8 blow, bup, betae, Ph, g0, umax, alpha, reparea, exfire

!! These variables are defined in params/fire
!
!    blow     : Lower biomass threshold (kgC m-2)
!    bup      : Upper biomass threshold (kgC m-2)
!    betae    : Extinction wetness content (dimensionless) 
!    Ph       : Fire occurrence probability due to human causes (fraction)
!    g0       : Spread fire rate dependence on wind for calm conditions (dimensionless)
!    alpha    : Factor of decay for spread fire rate dependence on wind (dimensionless)
!    umax     : Maximum fire spread rate (km h-1)
!    reparea  : Representative area of 1000 km2
!    exfire   : Fire extinguishing probability (fraction)

      real*8, dimension(:), allocatable :: Pbio, Pmoi, Pign, Pfire, srate, lbratio, abday, burnfrac, pfireyr, vegfrac

!    Pbio(npoi)           ! fire probability dependence on biomass (fraction)
!    Pmoi(npoi)           ! fire probability dependence on wetness (fraction)
!    Pign(npoi)           ! fire probability dependence on ignition (lightning and human causes) (fraction)
!    Pfire(npoi)          ! fire probability = Pbio*Pmoi*Pi (fraction)
!    srate                ! fire spread rate (km/h)
!    lbratio              ! length to breadth ratio of fire (dimensionless)
!    abday                ! burned area in 1 day (km2)
!    burnfrac             ! average burned fraction  (fraction)
!    pfireyr              ! annual mean fire probability (fraction)

      real*8, dimension(:,:), allocatable :: burnpft

!    burnpft              ! burned area for each PFT (fraction)

end module inland_comfire
