#include "inland_config.h"
module inland_comwork
      use inland_parameters
      implicit none
      public
      save

! -------
! comwork
! -------
!
! this file holds work space arrays and the longitude
! and latitude vectors.  The value of ndim3 should be nlon*nlat times
! the largest 3rd dimension in compar (max of nband,nsoilay,nsnolay,npft).
!   Ndim3 is now being defined on inland_alloc, so actual calculation of
! max(nband,nsoilay,nsnolay,npft) is done.

!    The variable ndim4 is calculated on main_offline subroutine. For the 
! offline model, the needed information is provided as compile-time macros on
! inland_parameters module. For now, this is how it works, but dimensions should be 
! allowed to change without the need of recompiling the code.
!    For that matter, ndim4 will be calculated on main_offline in order to fa-
! cilitate calculating it once values are provided by the inland.infile curren-
! tly read on main_offline subroutine.
      integer ndim2, ndim3, ndim4

      real*8, parameter :: ocean = 9.e+20
 
! use to store short names
      character aname*21

! lonindex(npoi): i index of nth point in land+sea array (i,j)
! latindex(npoi): j index of nth point in land+sea array (i,j)
! lonscale(nlon): longitude of nth point in degrees east [nlon]
! latscale(nlat): latitude of nth point in degrees north [nlat]
      integer, allocatable :: lonindex(:), latindex(:)     
      real*8, dimension(:), allocatable :: lonscale(:), latscale(:)

! work(ndim2): work space big enough for one full grid [ndim2]
! cdummy(ndim3): work space big enough for npft grids [ndim3]
! cdummyint(ndim3): work space big enough for npft grids [ndim3] (integer)
! buffer(npoi): work space big enough for all grid points [npoi]
! buffer1(npoi1): work space big enough for all non-tile grid points [npoi1]
      real*8, dimension(:), allocatable :: work(:)
      real*8, dimension(:), allocatable, target :: cdummy(:)
      integer, dimension(:), allocatable, target :: cdummyint(:)
      real*8, dimension(:), allocatable, target :: buffer(:)
      real*8, dimension(:), allocatable, target :: buffer1(:)

! seedvec(npoi): seeds used by random generator, unique for each point
      integer,  dimension(:), allocatable :: seedvec
 
end module inland_comwork
