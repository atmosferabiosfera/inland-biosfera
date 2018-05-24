#include "inland_config.h"
module inland_lsmmapib
      implicit none
      public
      save

! The land surface model works by gathering all the land points on the
! CCM3 [plon] x [plat] grid into a vector of [lpt] land points.
!
! [ixy] and [jxy] are indices for the mapping: [plon] x [plat] grid <->
! [lpt] vector of land points 

! Multiplicity of land points in arrays below
      integer, dimension(:), allocatable :: &
               ixy,   &! longitude index of lpt land points: 1 to plon
               jxy,   &! latitude index of lpt land points: 1 to plon
               begpnt,&! first point in "big" vector to process
               numpnt,&! number of points in "big" lpt vector to process
               procj, &! land point process number corres. to latitude decomp.
               procn   ! land point process number corres. to land decomp.
      real*8, dimension(:), allocatable :: &
              lati,   &! latitude of land point [radians]
              loni     ! longitude of land point [radians]

! Analogs on actual land mesh (ending in 1)
      integer, dimension(:), allocatable :: ixy1, jxy1, procj1, procn1
      real*8, dimension(:), allocatable :: lati1, loni1
end module inland_lsmmapib
