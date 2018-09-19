#include "inland_config.h"
! ----------------------------------------------------------------------------
subroutine initialcrop (irestart)
! ----------------------------------------------------------------------------
!
! common blocks
!
      use inland_parameters
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
!
      integer iy, j, ic
!
! reset variables after harvest date and data is written
! to io.f and crop diagnostic output files
! have to make sure there is no memory in the system
! of upper canopy or lower canopy vegetation when crops
! replace natural vegetation!
!
!      do j=scpft,ecpft
!      thrlai(:,j)=0.0
!      peaklai(:,j)=0.0

      cumlvs(:,:)=0.0
      dpgf(:,:)=0.0

!      totnuptake(:,j)=0.0
!      tnplant(:,j)=0.0
!      totnfix(:,j)=0.0
!      idpp(:,j)=0.0
!      fixn(:,j)=0.0
!      gddplant(:,j)=0.0
!      crmclim(:,j)=0.0
!      crmact(:,j)=0.0
!      crmplant(:,j)=0.0
!      grainday(:,j)=9999.0
!      gddtsoi(:,j)=0.0
!      fertinput(:,j)=0.0
!      pstart(:,j)=999
!      cdays(:)=0.0
!      cropy(:)=0.0
!      ik(:)=0.0
      if (irestart .eq. 0)then
         aylprod(:,:)=0.0
         ayrprod(:,:)=0.0
         ayabprod(:,:)=0.0
         aybprod(:,:)=0.0
         cbiol(:,:)=0.0
         cbios(:,:)=0.0
         cbior(:,:)=0.0
         cbiow(:,:)=0.0
         cbiog(:,:)=0.0
         cropout(:,:,:)=0.
         plaimx(:,:)=0.0
         plai(:,:)=0.0
         biomass(:,:)=0.0
         htmx(:,:)=0.0
         harvidx(:,:)=0.0
         hui(:,:)=0.0
         aleaf(:,:)=0.0
         aroot(:,:)=0.0
         astem(:,:)=0.0
         arepr(:)=0.0
         abunch(:,:)=0.0
         awood(:,:)=0.0
         leafout(:,:)=0.0
      endif
!     enddo

! grid cell variables - not dependent on pft
!
!      sai(:,j)=0.0
      fu(:)=0.0
!      lai(:,j)=0.0
!      zbot(:,j)=0.0
!      ztop(:,j)=0.0
!      totbiou(:)=0.0
!      totbiol(:)=0.0
!      totlaiu(:)=0.0
!      totlail(:)=0.0
!      vf(:)=0.0
!      ncyears(:)=0.0

!
      return
!
end subroutine initialcrop
