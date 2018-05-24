#include "inland_config.h"
! --------------------------------------------------------------------
subroutine cellbox(lat, latp1, lon, lonp1, edgen, edges, edgew, edgee, numlon, longxy, latixy, lats , lonw)
! --------------------------------------------------------------------

! ------------------------ code history ------------------------------
! source file       : cellbox.F
! purpose           : southern and western edges of grid cells
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
      real*8 edgen              !northern edge of grid (degrees)
      real*8 edges              !southern edge of grid (degrees)
      real*8 edgew              !western edge of grid (degrees)
      real*8 edgee              !eastern edge of grid (degrees)
      integer numlon(lat)     !number of grid cells per latitude strip
c     real*8 longxy(lon,lat)    !longitude at center of grid cell
c     real*8 latixy(lon,lat)    !latitude at center of grid cell
      real*8 longxy(lon)    !longitude at center of grid cell
      real*8 latixy(lat)    !latitude at center of grid cell

! ------------------- output variables ----------------------------
      real*8 lats(latp1)        !grid cell latitude, southern edge (degrees)
      real*8 lonw(lonp1,lat)    !grid cell longitude, western edge (degrees)

! ------------------- local variables -----------------------------
      integer i               !longitude index
      integer j               !latitude index
      real*8 dx                 !cell width

! --------------------------------------------------------------------
! Latitudes -- southern edges for each latitude strip. The southern
! and northern edges of latitude strip j are:
!        southern = lats(j  )
!        northern = lats(j+1)
! Hence, [lats] is dimensioned lats(lat+1)
! --------------------------------------------------------------------     
      lats(1) = edges 
      do j = 2, lat
c        lats(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
         lats(j) = (latixy(j-1) + latixy(j)) / 2.
      end do
      lats(lat+1) = edgen 

! --------------------------------------------------------------------
! Longitudes -- western edges. Longitudes for the western edge of the 
! cells must increase continuously and span 360 degrees. Three types of 
! grids are valid:
!
! 1: grid starts at Dateline  with western edge on Dateline
! 2: grid starts at Greenwich with western edge on Greenwich
! 3: grid starts at Greenwich with center on Greenwich
!
! For Grid 1 (Dateline)          , western edges range from:  -180 to 180
! For Grid 2 (Greenwich)         , western edges range from:     0 to 360
! For Grid 3 (Greenwich centered), western edges range from: -dx/2 to -dx/2 + 360 
!
! Western edges correspond to [longxy] (longitude at center of cell) for
! Grid 1 only. In this case, western edges range from -180 to 180 with
! negative longitudes west of Greenwich. Hence, this is the preferred
! grid type. Grids 2 and 3 are supported because some data sets start
! at Greenwich rather than Dateline (Grid 2) and the NCAR CCM starts
! at Greenwich, centered on Greenwich (Grid 3). 
!
! Partial grids that do not span 360 degrees are allowed so long as they
! have the convention of Grid 1 with 
!      western edge of grid: >= -180 and < 180
!      eastern edge of grid: > western edge  and <= 180
!
! [lonw] must be dimensioned lonw(lon+1,lat) because each latitude
! strip can have variable longitudinal resolution
! --------------------------------------------------------------------

      do j = 1, lat
         dx = (edgee - edgew) / numlon(j)
! western edge of first grid cell

c         lonw(1,j) = longxy(1,j) - dx/2.
         lonw(1,j) = longxy(j) - dx/2.

! remaining grid cells
         do i = 2, numlon(j)+1
            lonw(i,j) = lonw(1,j) + (i-1)*dx
         end do

! set unused longitudes to non-valid number
         do i = numlon(j)+2, lon
            lonw(i,j) = -999.
         end do
      end do
      return
end subroutine cellbox
