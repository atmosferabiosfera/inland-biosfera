#include "inland_config.h"
module inland_comtex
      implicit none
      public
      save
! ------
! comtex
! ------
!
! ----------------------------------------
! Soil texture-related parameters for inland
! ----------------------------------------

      real*8, dimension(:), allocatable ::  porosdat, sfielddat, swiltdat,  &
                                            bexdat, suctiondat, hydrauldat, &
                                            cpwfdat
!
!      porosdat(ndat),   ! porosity volume fraction
!      sfielddat(ndat),  ! field capacity volume fraction
!      swiltdat(ndat),   ! wilting point volume fraction
!      bexdat(ndat),     ! Campbell moisture-release b exponent
!      suctiondat(ndat), ! Air entry potential (m-H20)
!      hydrauldat(ndat), ! saturated hydraulic conductivity (m s-1)
!      cpwfdat(ndat)     ! Kaiyuan Li for Green-Ampt, capillary pressure at wetting front (m)  

      real*8, dimension(:,:), allocatable ::  texdat
!
!      texdat(3, ndat),  ! sand/silt/clay fractions
!
end module inland_comtex
