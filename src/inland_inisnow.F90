#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine inisnow 
! ---------------------------------------------------------------------
! does initialization for snow model
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsoi
      use inland_comsno
      use inland_comcrop, only:isimagro

      implicit none
     
     if(isimagro .gt. 0) then 
! ---------------------------------------------------------------------
! rhos is density of snow
      rhos = 0.15 * rhow
! consno is thermal conductivity of snow
      consno = 0.20

! hsnotop is "adaptive-grid" thickness of top snow layer
      hsnotop = 0.1

! hsnomin is minimum total snow thickness. total thickness
! is constrained to hsnomin for less than 100% cover. (hsnomin
! should be ge nsnolay*hsnotop for vadapt to work properly.)
      hsnomin = max (dble(0.3), nsnolay * hsnotop)

! fimin and fimax are minimum and maximum snowcover fractions
      fimin = 0.00002 * (dtime / 1800.) * (0.3 / hsnomin)
      fimax = 1.000
! z0sno is roughness lenth of snow cover
      z0sno = 0.0005


     else

! rhos is density of snow
      rhos = 0.25 * rhow

! consno is thermal conductivity of snow
      consno = 0.20

! hsnotop is "adaptive-grid" thickness of top snow layer
      hsnotop = 0.05

! hsnomin is minimum total snow thickness. total thickness
! is constrained to hsnomin for less than 100% cover. (hsnomin
! should be ge nsnolay*hsnotop for vadapt to work properly.)
      hsnomin = max (dble(0.15), nsnolay * hsnotop)

! fimin and fimax are minimum and maximum snowcover fractions
      fimin = 0.00002 * (dtime / 1800.) * (0.15 / hsnomin)
      fimax = 1.000

! z0sno is roughness lenth of snow cover
      z0sno = 0.0005

     endif


      return
end subroutine inisnow
