#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine existence
! ---------------------------------------------------------------------
!
! this routine determines which plant functional types (pft's) are allowed
! to exist in each gridcell, based on a simple set of climatic criteria
!
! the logic here is based on the biome3 model of haxeltine and prentice
!
! plant functional types:
!
! 1)  tropical broadleaf evergreen trees
! 2)  tropical broadleaf drought-deciduous trees
! 3)  warm-temperate broadleaf evergreen trees
! 4)  temperate conifer evergreen trees
! 5)  temperate broadleaf cold-deciduous trees
! 6)  boreal conifer evergreen trees
! 7)  boreal broadleaf cold-deciduous trees
! 8)  boreal conifer cold-deciduous trees
! 9)  evergreen shrubs
! 10) deciduous shrubs
! 11) warm (c4) grasses
! 12) cool (c3) grasses
! 13) c3 crop (soybean)
! 14) c4 crop (corn)
! 15) c3 crop (wheat)
! 16) c4 crop (sugarcane)


      use inland_parameters
      use inland_control, only: overveg
      use inland_comcrop
      use inland_compft
      use inland_comveg
      use inland_comsum
      use inland_comnitr
      use inland_combcs, only: xinveg

      implicit none
!-----------------------------------------------------------------------
!
! local variables
!
      integer i,j, & ! loop counters
              inveg  ! vegetation type
!
      do 100 i = lbeg, lend

! determine which plant types can exist in a given gridcell

        exist(i,1)  = 0.
        exist(i,2)  = 0.
        exist(i,3)  = 0.
        exist(i,4)  = 0.
        exist(i,5)  = 0.
        exist(i,6)  = 0.
        exist(i,7)  = 0.
        exist(i,8)  = 0.
        exist(i,9)  = 0.
        exist(i,10) = 0.
        exist(i,11) = 0.
        exist(i,12) = 0.
        exist(i,13) = 0.
        exist(i,14) = 0.
        exist(i,15) = 0.
        exist(i,16) = 0.
!
!*** DTP 2001/06/07: Modified version of above code reads in PFT
! existence criteria from external parameter file "params.veg"
!    These are copied here for reference.... 
!------------------------------------------------------------------
!  TminL    TminU    Twarm    GDD    PFT
!------------------------------------------------------------------
!    0.0   9999.0   9999.0   9999  !   1
!    0.0   9999.0   9999.0   9999  !   2
!  -10.0      0.0   9999.0   9999  !   3
!  -45.0      0.0   9999.0   1200  !   4
!  -45.0      0.0   9999.0   1200  !   5
!  -57.5    -45.0   9999.0    350  !   6
!  -57.5    -45.0   9999.0    350  !   7
! 9999.0    -45.0   9999.0    350  !   8
! 9999.0   9999.0   9999.0    100  !   9
! 9999.0   9999.0   9999.0    100  !  10
! 9999.0   9999.0     22.0    100  !  11
! 9999.0   9999.0   9999.0    100  !  12
!------------------------------------------------------------------
!
! 1) tropical broadleaf evergreen trees
!
!  - tcmin > 0.0
         if (tcmin(i).gt.TminL(1)) exist(i,1) = 1.0

! 2) tropical broadleaf drought-deciduous trees
!
!  - tcmin > 0.0
         if (tcmin(i).gt.TminL(2)) exist(i,2) = 1.0

! 3) warm-temperate broadleaf evergreen trees
!
!  - tcmin <   0.0 and
!  - tcmin > -10.0
         if ((tcmin(i).lt.TminU(3)).and.(tcmin(i).gt.TminL(3))) exist(i,3) = 1.0

! 4) temperate conifer evergreen trees
!
!  - tcmin <    0.0 and
!  - tcmin >  -45.0 and
!  - gdd5  > 1200.0
         if ((tcmin(i).lt.TminU(4)).and.(tcmin(i).gt.TminL(4)).and. &
             (gdd5(i).gt.GDD(4))) exist(i,4) = 1.0

! 5) temperate broadleaf cold-deciduous trees
!
!  - tcmin <    0.0 and
!  - tcmin >  -45.0 and
!  - gdd5  > 1200.0
         if ((tcmin(i).lt.TminU(5)).and.(tcmin(i).gt.TminL(5)).and. &
             (gdd5(i).gt.GDD(5))) exist(i,5) = 1.0

! 6) boreal conifer evergreen trees
!
!  - tcmin <  -45.0 or gdd5 < 1200.0, and
!  - tcmin >  -57.5 and
!  - gdd5  >  350.0
         if (((tcmin(i).lt.TminU(6)).or.(gdd5(i).lt.GDD(4))).and. &
              (tcmin(i).gt.TminL(6)).and.(gdd5(i).gt.GDD(6))) exist(i,6) = 1.0

! 7) boreal broadleaf cold-deciduous trees
!
!  - tcmin <  -45.0 or gdd5 < 1200.0, and
!  - tcmin >  -57.5 and
!  - gdd5  >  350.0
         if (((tcmin(i).lt.TminU(7)).or.(gdd5(i).lt.GDD(5))).and. &
             (tcmin(i).gt.TminL(7)).and.(gdd5(i).gt.GDD(7))) exist(i,7) = 1.0

! 8) boreal conifer cold-deciduous trees
!
!  - tcmin <  -45.0 or gdd5 < 1200.0, and
!  - gdd5  >  350.0
         if (((tcmin(i).lt.TminU(8)).or.(gdd5(i).lt.GDD(4))).and. &
              (gdd5(i).gt.GDD(8))) exist(i,8) = 1.0

! 9) evergreen shrubs
!
!  - gdd0 > 100.0
         if (gdd0(i).gt.GDD(9)) exist(i,9) = 1.0

! 10) deciduous shrubs
!
!  - gdd0 > 100.0
         if (gdd0(i).gt.GDD(10)) exist(i,10) = 1.0

! 11) warm (c4) grasses
!
!  - tw   >  22.0 and
!  - gdd0 > 100.0
         if ((tw(i).gt.Twarm(11)).and.(gdd0(i).gt.GDD(11))) exist(i,11) = 1.0

! 12) cool (c3) grasses
!
!  - gdd0 > 100.0
         if (gdd0(i).gt.GDD(12)) exist(i,12) = 1.0

!
!
! == C. Kucharik 6.12.01 ==
!
! if override natural vegetation competition (overveg = 1)
! this code is used to override existence parameterizations for potential
! vegetation distribution based on climatic constraints. Instead
! we only allow PFTs to compete in each grid cell
! based on land cover dataset and classification found in that region 
! override those pfts that are not desired but might have exist = 1.0
! from above initialization - this essentially limits vegetation competition
! during spin-up periods so that vegetation growing there is confined to
! what is typically observed today (potential vegetation).  If doing 
! climate change scenarios, overveg should be set to 0 so full 
! vegetation dynamics are used, if desired. 
!
         if(isimagro .gt. 0)then
            overveg=1
         else
            overveg=0
         endif

         if (overveg .eq. 1) then
!
            inveg = nint(xinveg(i))

!          
! wherever vegetation is allowed to exist, check the rules from above 
!
! tropical deciduous
!
            if (inveg .eq. 2) then
               exist(i,2) = 1.0
!
               exist(i,1) = 0.0
               do 70 j = 3,16
                  exist(i,j) = 0.0
 70            continue

!
! temperate conifers
!
            else if (inveg .eq. 4) then
               exist(i,4) = 1.0
!
               do 71 j = 1,3
                  exist(i,j) = 0.0
 71            continue
!
               do 72 j = 5,16
                  exist(i,j) = 0.0
 72            continue
!
! temperate deciduous
!
            else if (inveg .eq. 5) then
               exist(i,5) = 1.0
!
               do 73 j = 1,4
                  exist(i,j) = 0.0
 73            continue
!
               do 74 j = 6,16
                  exist(i,j) = 0.0
 74            continue
!
! boreal conifers
!
            else if (inveg .eq. 6) then
               exist(i,6) = 1.0
!
               do 75 j = 1,5
                  exist(i,j) = 0.0
 75            continue
!
               do 76 j = 7,16
                  exist(i,j) = 0.0
 76            continue
!
! boreal deciduous 
!
            else if (inveg .eq. 7) then
               do 77 j = 7, 8
                  exist(i,j) = 1.0
 77            continue
!
               do 78 j = 1,6
                  exist(i,j) = 0.0
 78            continue
!
               do 79 j = 9,16
                  exist(i,j) = 0.0
 79            continue
!
! mixed forest exception:
!
! let existence rules determine whether temperate or boreal species
! can exist at this location
! 
            else if (inveg .eq. 8) then
!              do 69 j = 4, 8
!                 exist(i,j) = 1.0
! 69           continue
!
               do 81 j = 1,3
                  exist(i,j) = 0.0
 81            continue
!
               do 82 j = 9,16
                  exist(i,j) = 0.0
 82            continue
!
! savanna
!
            else if (inveg .eq. 9) then
               exist(i,5) = 1.0
               exist(i,11) = 1.0
               exist(i,12) = 1.0
!
               do 83 j = 1, 4
                  exist(i,j) = 0.0
 83            continue
!
               do 84 j = 6, 10 
                  exist(i,j) = 0.0
 84            continue
!
               do 85 j = 13,16
                  exist(i,j) = 0.0
 85            continue
!
! grassland
!
            else if (inveg .eq. 10) then
               exist(i,11) = 1.0
               exist(i,12) = 1.0
!
               do 86 j = 1, 10 
                  exist(i,j) = 0.0
 86            continue
!
               do 87 j = 13, 16 
                  exist(i,j) = 0.0
 87            continue
!
! dense shrubland 
!
            else if (inveg .eq. 11) then
               exist(i,9)  = 1.0
               exist(i,10) = 1.0
               exist(i,11) = 1.0
               exist(i,12) = 1.0
!
               do 88 j = 1, 8 
                  exist(i,j) = 0.0
 88            continue
!
               do 89 j = 13, 16 
                  exist(i,j) = 0.0
 89            continue
!
! open shrubland
!
            else if (inveg .eq. 12) then
!
               do 90 j = 9, 12 
                  exist(i,j) = 1.0
 90            continue
!
               do 91 j = 1, 8 
                  exist(i,j) = 0.0
 91            continue
!
               do 92 j = 13, 16 
                  exist(i,j) = 0.0
 92            continue
!
! tundra 
!
            else if (inveg .eq. 13) then
!
               do 93 j = 9, 12 
                  exist(i,j) = 1.0
 93            continue
!
               do 94 j = 1, 8
                  exist(i,j) = 0.0
 94            continue
!
               do 95 j = 13, 16 
                  exist(i,j) = 0.0
 95            continue
            endif
         endif ! override original natural vegetation competition
!
! == C. Kucharik ==
! cropping systems - are currently planted everywhere


! if landusetype is "natural" (1) then only natural vegetation rules apply,
!   but if landusetype >= 2 then some pfts may not be allowed to exist

! if a point has landusetype "cropland" (2) and a specific croptype, then only the crop pft can exist
! else (because croptype is not known or landusetype is pasture/urban), use values in landusepftmap 
!   (which is built from the last table in params/vegetation)

         if ( landusetype(i) .gt. 1 ) then
            if ( (landusetype(i) .eq. 2) .and. (croptype(i) .ge. scpft) ) then
               ! only crop pft can exist if a crop is defined at point
               exist(i,:) = 0.
               exist(i,nint(croptype(i))) = 1.
            else
               ! generic case, some pfts may not exist
               do j = 1, npft
                  if ( landusepftmap(nint(landusetype(i)),j) .eq. 0. ) exist(i,j) = 0.
               end do
            end  if
         end if

100   continue

      return
end subroutine existence
