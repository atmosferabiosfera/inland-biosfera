#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine cropupdate(jday)
! ---------------------------------------------------------------------
!
! common blocks
!
      use inland_parameters
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
      use inland_comsno
      use inland_compft
!
!
!     real*8 laimx(npft)
!
      real*8 avglail
!	     xminlai,   &
!
      integer j, i, jday
!
! begin global grid

      do 100 i = lbeg, lend
         if (croptype(i) .ge. scpft) then
            do 90 j = scpft, ecpft 
               if (exist(i,j) .eq. 1.0 .and. croplive(i,j) .eq. 1.0) then
                  cbiol(i,j) = max((exist(i,j) * xminlai/specla(i,j)),cbiol(i,j))
                  if (plai(i,j) .gt. plaimx(i,j)) plaimx(i,j) = plai(i,j)
               endif
90          continue

! crop canopy single sided leaf area index (area-weighted)

            avglail = plai(i,13) + plai(i,14) + plai(i,15) + plai(i,16)           

! crop canopy fractions

            frac(i,13) = plai(i,13) / max (avglail, epsilon)
            frac(i,14) = plai(i,14) / max (avglail, epsilon)
            frac(i,15) = plai(i,15) / max (avglail, epsilon)
            frac(i,16) = plai(i,16) / max (avglail, epsilon)

! calculate total crop leaf are index

            totlail(i) = plai(i,13) + plai(i,14) + plai(i,15) + plai(i,16)
            fl(i) = totlail(i) / 1.0
            fl(i) = max(0.025, min(0.975, fl(i)))

! calculate the crop canopy leaf area index using the fractional vegetation cover

            lai(i,1) = avglail / fl(i)

! C. Kucharik  04.02.01

            greenfracl(i) = 0.0
            do 80 j = scpft, ecpft
               greenfracl(i) = greenfracl(i) + frac(i,j) * grnfraccrop(i,j)
80          continue
!
! calculate total crop canopy biomass
!
            totbiol(i) = biomass(i,13) + biomass(i,14) + &
			 biomass(i,15) + biomass(i,16)
            zbot(i,1) = 0.02
            ztop(i,1) = plai(i,13) / lai(i,1) * ztopmxsoy * &
                        (min(plai(i,13) / (laimx(13)-1.0),1.0))**2 + &
                        plai(i,14) / lai(i,1) * ztopmxmze * &
                        (min(plai(i,14) / (laimx(14)-1.5),1.0))**2 + &
                        plai(i,15) / lai(i,1) * ztopmxwht * &
                        (min(plai(i,15) / (laimx(15)-1.0),1.0))**2 + &
                        plai(i,16) / lai(i,1) * ztopmxsgc * &
			            min((rm(i) / 50.0), 1.0) * &
                        (min(plaimx(i,16) / (laimx(16)),1.0))**2 
            htmx(i,1) = max(htmx(i,1),ztop(i,1))
            ztop(i,1) = max(0.05, max(htmx(i,1),ztop(i,1)))

! calculate stem area index for crops

            sai(i,1) = 0.20 * plai(i,13) + 0.10 * plai(i,14) + 0.20 * &
		       plai(i,15) + 0.10 * plai(i,16) 

! calculate annual aboveground npp total

            ayanpptot(i) = ayanpp(i,13) + ayanpp(i,14) + &
			   ayanpp(i,15) + ayanpp(i,16)

! end of loop

         endif  ! existence
100   continue

!! return to main program
 
      return
!
end subroutine cropupdate
