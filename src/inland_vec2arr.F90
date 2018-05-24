#include "inland_config.h"

! ---------------------------------------------------------------------
subroutine vec2arr(vect, array, ilpt, nodata)
! puts vector of land-only points back into array
! ilpt argument is the subgrid #
! for original version without tiles, use vec2arr(vect, array, 1, ocean) instead
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comwork, only: lonindex, latindex
      use inland_subgrid

      implicit none
!------------------------------Arguments--------------------------------
      real*8, intent(out) :: array(plona,plata)
      real*8, intent(in)  :: vect(lbeg:lend)
      real*8, intent(in)  :: nodata
      integer, intent(in) :: ilpt

      integer i,k !i is in array, k is in vect (counting ilpt)

      array(:,:)=nodata
      do i=lbeg,lend/mlpt
         k=subgrid_get_index(i,ilpt)
         if ( k .eq. 0 ) then
            write (*,*) 'ERROR in vec2arr2(), got subgrid index',k,'with i=',i,' ilpt=',ilpt
         end if
         array(lonindex(i),latindex(i)) = vect(k)
         !print *,i,lonindex(i),latindex(i),vect(k)
      end do

      return
end subroutine vec2arr

! ---------------------------------------------------------------------
subroutine vec2arr2d(vect, array, ilpt, nodata, ixdim, nxdim)
! version that takes a 2D array as input, used in writevar
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comwork, only: lonindex, latindex
      use inland_subgrid

      implicit none
!------------------------------Arguments--------------------------------
      real*8, intent(out) :: array(plona,plata)
      real*8, intent(in)  :: vect(lbeg:lend,nxdim)
      real*8, intent(in)  :: nodata
      integer, intent(in) :: ilpt,ixdim,nxdim

      integer i,j,k !i is in array, k is in vect (counting ilpt)

      array(:,:)=nodata
      do i=lbeg,lend/mlpt
         k=subgrid_get_index(i,ilpt)
         if ( k .eq. 0 ) then
            write (*,*) 'ERROR in vec2arr2(), got subgrid index',k,'with i=',i,' ilpt=',ilpt
         end if
         array(lonindex(i),latindex(i)) = vect(k,ixdim)
         !print *,i,lonindex(i),latindex(i),vect(k)
      end do

      return
end subroutine vec2arr2d