#include "inland_config.h"
!-------------------------------------------------------------------------
integer function textcls(msand,mclay)
! adapted by cjk 01/11/01
!-------------------------------------------------------------------------
! |
! |                         T R I A N G L E
! | Main program that calls WHAT_TEXTURE, a function that classifies soil
! | in the USDA textural triangle using sand and clay %
! +-----------------------------------------------------------------------
! | Created by: aris gerakis, apr. 98 with help from brian baer
! | Modified by: aris gerakis, july 99: now all borderline cases are valid
! | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
! |              main program
! +-----------------------------------------------------------------------
! | COMMENTS
! | o Supply a data file with two columns, in free format:  1st column sand,
! |   2nd column clay %, no header.  The output is a file with the classes.
! +-----------------------------------------------------------------------
! | You may use, distribute and modify this code provided you maintain
! ! this header and give appropriate credit.
! +-----------------------------------------------------------------------
!
! input variables
      integer  msand,&
               mclay

! local variables
      logical inpoly

      real*8  silty_loam(1:7,1:2),      &
              sandy(1:7,1:2),           &
              silty_clay_loam(1:7,1:2), &
              loam(1:7,1:2),            &
              clay_loam(1:7,1:2),       &
              sandy_loam(1:7,1:2),      &
              silty_clay(1:7,1:2),      &
              sandy_clay_loam(1:7,1:2), &
              loamy_sand(1:7,1:2),      &
              clayey(1:7,1:2),          &
              sandy_clay(1:7,1:2)

! ---------------------------------------------------------------------
! initalize polygon coordinates:
! each textural class reads in the sand coordinates (1,7) first, and
! then the corresponding clay coordinates (1,7)

!     data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
!
! because we do not have a separate silt category, have to redefine the
! polygon boundaries for the silt loam  
      data sandy           /85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
      data loamy_sand      /70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
      data sandy_loam      /50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 15, 0/
      data loam            /43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
      data silty_loam      /0, 0, 23, 50, 0, 0, 0, 0, 27, 27, 0, 0, 0, 0/ 
      data sandy_clay_loam /52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35, 20, 0, 0/
      data clay_loam       /20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data silty_clay_loam /0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data sandy_clay      /45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/
      data silty_clay      /0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
      data clayey          /20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/

! polygon coordinates  
!
!     sand
!
!     >  85, 90, 100, 0, 0, 0, 0,       ! sand
!     >  70, 85, 90, 85, 0, 0, 0,       ! loamy sand
!     >  50, 43, 52, 52, 80, 85, 70,    ! sandy loam
!     >  43, 23, 45, 52, 52, 0, 0,      ! loam
!     >   0, 0, 23, 50, 0, 0, 0,        ! silt loam (combined with silt)
!     >  52, 45, 45, 65, 80, 0, 0,      ! sandy clay loam
!     >  20, 20, 45, 45, 0, 0, 0,       ! clay loam
!     >   0, 0, 20, 20, 0, 0, 0,        ! silty clay loam
!     >  45, 45, 65, 0, 0, 0, 0,        ! sandy clay
!     >   0, 0, 20, 0, 0, 0, 0,         ! silty clay 
!     >  20, 0, 0, 45, 45, 0, 0         ! clay
!
!      clay
!
!     > 0, 10, 0, 0, 0, 0, 0,           ! sand
!     > 0, 15, 10, 0, 0, 0, 0,          ! loamy sand
!     > 0, 7, 7, 20, 20, 15, 0,         ! sandy loam 
!     > 7, 27, 27, 20, 7, 0, 0,         ! loam
!     > 0, 27, 27, 0, 0, 0, 0,          ! silt loam (combined with silt)
!     > 20, 27, 35, 35, 20, 0, 0,       ! sandy clay loam
!     > 27, 40, 40, 27, 0, 0, 0,        ! clay loam
!     > 27, 40, 40, 27, 0, 0, 0,        ! silty clay loam
!     > 35, 55, 35, 0, 0, 0, 0,         ! sandy clay
!     > 40, 60, 40, 0, 0, 0, 0,         ! silty clay
!     > 40, 60, 100, 55, 40, 0, 0       ! clay
!
! +-----------------------------------------------------------------------
! | figure out what texture grid cell and layer are part of  
! | classify a soil in the triangle based on sand and clay %
! +-----------------------------------------------------------------------
! | Created by: aris gerakis, apr. 98
! | Modified by: aris gerakis, june 99.  Now check all polygons instead of
! | stopping when a right solution is found.  This to cover all borderline 
! | cases.
! +-----------------------------------------------------------------------
!
! find polygon(s) where the point is.  
      textcls = 0 
      if (msand .gt. 0.0 .and. mclay .gt. 0.0) then
         if (inpoly(sandy, 3, msand, mclay)) then
            textcls = 1      ! sand
         endif
         if (inpoly(loamy_sand, 4, msand, mclay)) then
            textcls = 2      ! loamy sand
         endif
         if (inpoly(sandy_loam, 7, msand, mclay)) then
            textcls = 3      ! sandy loam
         endif
         if (inpoly(loam, 5, msand, mclay)) then
            textcls = 4      ! loam
         endif
         if (inpoly(silty_loam, 4, msand, mclay)) then
            textcls = 5      ! silt loam
         endif
         if (inpoly(sandy_clay_loam, 5, msand, mclay)) then
            textcls = 6      ! sandy clay loam
         endif
         if (inpoly(clay_loam, 4, msand, mclay)) then
            textcls = 7      ! clay loam
         endif
         if (inpoly(silty_clay_loam, 4, msand, mclay)) then
            textcls = 8      ! silty clay loam
         endif
         if (inpoly(sandy_clay, 3, msand, mclay)) then
            textcls = 9      ! sandy clay
         endif
         if (inpoly(silty_clay, 3, msand, mclay)) then
            textcls = 10     ! silty clay
         endif
         if (inpoly(clayey, 5, msand, mclay)) then
            textcls = 11     ! clay
         endif
      endif

      if (textcls .eq. 0) then
         textcls = 5         ! silt loam
!        write (*, 1000) msand, mclay
!1000    format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
      endif
end function textcls
