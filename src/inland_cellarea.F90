#include "inland_config.h"
! --------------------------------------------------------------------
subroutine cellarea (lat, latp1, lon, lonp1, edgen, edges, &
                     edgew, edgee, numlon, lats, lonw, area)
! --------------------------------------------------------------------

! ------------------------ code history ------------------------------
! source file       : cellarea.F
! purpose           : area of grid cells (square kilometers)
! date first created: March 1996 - lsm version 1 (dataset generation code)
! by whom           : Gordon Bonan
! date last revised : December 1998 - lsm version 2
! by whom           : Gordon Bonan
! --------------------------------------------------------------------

      implicit none

! ------------------------ input variables ---------------------------
      integer lat             !dimension: number of latitude points
      integer latp1           !dimension: lat + 1
      integer lon             !dimension: number of longitude points
      integer lonp1           !dimension: lon + 1
      real*8 edgen            !northern edge of grid (degrees)
      real*8 edges            !southern edge of grid (degrees)
      real*8 edgew            !western edge of grid (degrees)
      real*8 edgee            !eastern edge of grid (degrees)
      integer numlon(lat)     !number of grid cells per latitude strip
      real*8 lats(latp1)      !grid cell latitude, southern edge (degrees)
      real*8 lonw(lonp1,lat)  !grid cell longitude, western edge (degrees)
! --------------------------------------------------------------------

! ------------------------ output variables --------------------------
      real*8 area(lon,lat)    !cell area (km**2)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
      real*8 deg2rad          !pi/180
      real*8 re               !radius of earth (km)
      real*8 global           !summed area
      integer j               !latitude index
      integer i               !longitude index
      real*8 dx               !cell width: E-W
      real*8 dy               !cell width: N-S
      real*8 error            !true area for error check
! --------------------------------------------------------------------

      deg2rad = (4.*atan(1.)) / 180.
      re = 6371.227709
      global = 0.

      do j = 1, lat
         do i = 1, numlon(j)
            dx = (lonw(i+1,j) - lonw(i,j)) * deg2rad
            dy = sin(lats(j+1)*deg2rad) - sin(lats(j)*deg2rad)
            area(i,j) = dx*dy*re*re
            global = global + area(i,j)
         end do
      end do

! make sure total area from grid cells is same as area of grid
! as defined by its edges

      dx = (edgee - edgew) * deg2rad
      dy = sin(edgen*deg2rad) - sin(edges*deg2rad)
      error =  dx*dy*re*re
!
      if (abs(global-error)/error .gt. 0.0001) then
         write (6,*) 'CELLAREA error: correct area is ', error, &
                     ' but summed area of grid cells is ', global
         call endrun
      end if

      return
end subroutine cellarea

