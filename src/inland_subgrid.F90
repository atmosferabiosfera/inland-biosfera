#include "inland_config.h"
module inland_subgrid
!

      use inland_comwork ! only: all!

      implicit none
      public
      save

! ------
! subgrid
! ------
!
! tilefrac: the fraction of every subgrid pixel within the "parent" pixel - 1:npoi1
! waterfrac: the fraction of water cover in every "parent" pixel - 1:npoi
! tilefrac1: initial tilefracs, read from hr dataset(before redistributing out-of-range) - not allocated by default, unless hrmapfile is given
! tilefrac2: initial tilefracs (out-of-range, including water)
! itilechild: the indexes of all non-dominant subgrid pixels, has dims (npoi1,(maxsubgrid)) and value 0 if non-existent
! ntilechild: the number of subgrid tiles in this pixel - can also get it by finding children which tilefrac>0
! itileparent: the index of a subgrid "parent" grid
!              if tile is "dominant" then value is its index, or value 0 if non-existent - unused for now, set in subgrid_set_index
! hrmapnlon, hrmapnlat: dimensions of high res map
! hrmapres: resolution (x==y) of high res map
! ihtrileparent: the index of the tile (1:npoi1) for each high res map pixel
!   integer*2 has range from -32768 to 32767
!   use offset of -32768 to support global run (1:58920) - this is stored in ihrtileparentoffset
!   this takes around 300MB for sam run

      ! model-resolution variables
      integer maxsubgrid !ET TODO - parameter?
      real*8, dimension(:), allocatable :: tilefrac
      real*8, dimension(:), allocatable :: waterfrac
      integer, dimension(:,:), allocatable :: itilechild
      integer, dimension(:), allocatable :: ntilechild
      integer, dimension(:), allocatable :: itileparent

      ! hrmap variables, initialized in rdhrvegmap
      integer :: hrmapnlon, hrmapnlat
      real*8 hrmapres  
      real*8, dimension(:), allocatable :: hrmaplonvalues, hrmaplatvalues
      integer*2, dimension(:,:), allocatable :: ihrtileparent
      integer ihrtileparentoffset !32768
      real*8, dimension(:), allocatable :: tilefrac1, tilefrac2 !used in rdhrvegmap, tmp

contains

! ---------------------------------------------------------------------
integer function subgrid_get_index(inpoi,ilpt)
! ---------------------------------------------------------------------
! return the index of the subgrid tile #ilpt corresponding to dominant tile at i

      implicit none
      integer, intent(in) :: inpoi, ilpt

      subgrid_get_index = 0 ! default 0, non existent

      if ( ilpt .gt. 0 .and. ilpt .le. mlpt ) then
         subgrid_get_index = itilechild(inpoi,ilpt)
      end if

end function subgrid_get_index

! ---------------------------------------------------------------------
integer function subgrid_set_index(inpoi,ilpt,index)
! ---------------------------------------------------------------------
! set the index of the subgrid tile #ilpt corresponding to dominant tile at i

      implicit none
      integer, intent(in) :: inpoi, ilpt, index

      subgrid_set_index = 0 ! default 0, non existent

      subgrid_set_index = index
      itilechild(inpoi,ilpt) = index
      itileparent(index) = inpoi

end function subgrid_set_index

! ---------------------------------------------------------------------
subroutine subgrid_calculate_area(invar,outvar)
! ---------------------------------------------------------------------
! given a fractional variable (e.g. burnfrac), 
! calculate area value (e.g. burnarea) for each subgrid point

      implicit none
      real*8, dimension(:), intent(in) :: invar
      real*8, dimension(:), intent(out) :: outvar
      integer i,k,ilpt
      
      do i=1,npoi1
         do ilpt = 1, mlpt
            k=subgrid_get_index(i,ilpt)
            outvar(k) = invar(k) * garea(i) * tilefrac(k)
         end do
      end do

end subroutine subgrid_calculate_area


end module inland_subgrid
