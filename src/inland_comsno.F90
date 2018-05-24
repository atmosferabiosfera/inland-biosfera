#include "inland_config.h"
module inland_comsno

      implicit none
      public
      save
!
! ------
! comsno
! ------
!
      real*8 z0sno, rhos, consno, hsnotop, hsnomin, fimin, fimax
!
!     z0sno                ! roughness length of snow surface (m)
!     rhos                 ! density of snow (kg m-3)
!     consno               ! thermal conductivity of snow (W m-1 K-1)
!     hsnotop              ! thickness of top snow layer (m)
!     hsnomin              ! minimum total thickness of snow (m)
!     fimin                ! minimum fractional snow cover
!     fimax                ! maximum fractional snow cover
!
      real*8, dimension(:), allocatable :: fi
!
!     fi(npoi)             ! fractional snow cover
!
      real*8, dimension(:,:), allocatable :: tsno, hsno
!
!     tsno(npoi,nsnolay),  ! temperature of snow layers (K)
!     hsno(npoi,nsnolay)   ! thickness of snow layers (m)
!
end module inland_comsno
