#include "inland_config.h"
#include "inland_compar.h"

! ---------------------------------------------------------------------
subroutine arr2vec(array, vect)
! extracts land points from array and puts them into a land-only vector
! new version with tiles
! if mlpt > 1 it will replicate the values to all subgrid points
! ---------------------------------------------------------------------
      use inland_parameters

      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      real*8, intent(in)  :: array(plona,plata)
      real*8, intent(out) :: vect(lbeg:lend)

! Local variables
      integer ilpt
      
      ! replicate to all subgrid points - this is less efficient but code is shorter
      do ilpt = 1,mlpt
         call arr2vec_tile(array, vect, ilpt)
      end do

      return

end subroutine arr2vec

! ---------------------------------------------------------------------
subroutine arr2vec_tile(array, vect, ilpt)
! extracts land points from array and puts them into a land-only vector
! new version with tiles
! for original version without tiles, call arr2vec_tile(array, vect, 1)
! ilpt argument is the tile # to extract
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_combcs, only: lmask
      use inland_subgrid

      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      real*8, intent(in)  :: array(plona,plata)
      real*8, intent(out) :: vect(lbeg:lend)
      integer, intent(in) :: ilpt

! Local variables
      integer npts,i,j,k
      
      npts=0
      do 10 j= 1,plata
         do 20 i= 1,plona
            if (lmask(i,j) .eq. 1) then
               npts = npts + 1
               k = subgrid_get_index(npts,ilpt)
               if ( k .eq. 0 ) then
                  write (*,*) 'ERROR in arr2vec_tile(), got subgrid index',k,&
                       'with npts=',npts,' ilpt=',ilpt
               end if
               vect(k) = array(i,j)
            endif
20       continue
10    continue

      if (npts .ne. npoi/mlpt) then
         write (*,*) '*** ERROR: arr2vec: number of land points in lmask '// &
                     'does not match number of land points specified in '// &
                     'input parameters file.'
         stop 1
      endif
      return
end subroutine arr2vec_tile

