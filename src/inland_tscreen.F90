#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine tscreen (tscr,zscr,za,z1,z12,z2,z3,z34,z4,dispu,displ,ta,t12,t34,npt)
! ---------------------------------------------------------------------
!        Interpolates diagnostic screen-height temperature tscr at
!        height zscr.
!
! ---------------------------------------------------------------
     use inland_parameters
     implicit none
! -----------------------------------------------------------------------
!
! input variables
!
     integer :: npt             ! number of points in little vector
     real*8 :: zscr             ! refernce height
     real*8 ::  tscr(mpt),za(mpt),z1(mpt),z12(mpt),z2(mpt),z3(mpt),z34(mpt), &
                z4(mpt),dispu(mpt),displ(mpt),ta(mpt),t12(mpt),t34(mpt)
     ! variables displ, z2, z3 and z4 are unused

! local variables
     integer :: i
     real*8 :: w

     do 100 i = 1, npt
        if (zscr.gt.z1(i)) then
!         above upper canopy:
           w = log((zscr  - dispu(i)) / (z1(i) - dispu(i))) / &
               log((za(i) - dispu(i)) / (z1(i) - dispu(i)))
           tscr(i) = w * ta(i) + (1.-w) * t12(i)
        elseif (zscr.gt.z12(i)) then

!          within top half of upper canopy:
           tscr(i) = t12(i)
        else if (zscr.gt.z34(i)) then

!          between mid-points of canopies:
           tscr(i) = ( (zscr - z34(i)) * t12(i) + &
                       (z12(i) - zscr) * t34(i) ) &
                     / (z12(i)-z34(i))
        else

!          within or below lower canopy:
           tscr(i) = t34(i)
        endif
100   continue
      return
end subroutine tscreen
