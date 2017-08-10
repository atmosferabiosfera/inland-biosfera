#include "inland_config.h"
#include "inland_compar.h"
! ----------------------------------------------------------------------------
subroutine planting(irestart,iyrrestart,jday,ffact)
! ----------------------------------------------------------------------------
!
! subroutine to determine planting dates for crops based on observations
! and/or climatic constraints

! common blocks
      use inland_combcs
      use inland_parameters
      use inland_control,only:iday,imonth,iyear,iyear0
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
      use inland_comsno
!
! local variables
      real*8 pmmax(npft),      &
             pdmax(npft),      &
             sumhy(lbeg:lend,npft), &
             sumdp(lbeg:lend,npft), &
             yc, ffact

      integer jday,     &
              irestart, &
              iyrrestart, &
              i,j,k,n,  &
              npts,     &
              iy

      character*100 fertdata, &
                    header

      fertdata = 'hist.fert.1945.1996'   ! data file that contains historical
                                         ! changes in fertilizer usage across upper Midwest US  

! read in historical fertilizer data one time into array
! file has quantities in kg/ha, so these are changed to
! kg/m2 
!
! ** read in restart scenario also **
!      if ((iyear .eq. 1945 .or. iyear .eq. iyrrestart) .and. imonth .eq. 1 .and. iday .eq. 1) then
!         npts = 55
!         open(15, file = fertdata) 
!         read(15, *) header
!         do n = 1, npts 
!            read(15, *, end=20) iy, cfertmaize(iy),cfertsoy(iy)
!            cfertsgc(iy)=cfertmaize(iy)
!            cfertmaize(iy) = cfertmaize(iy) * 1e-04
!            cfertsgc(iy)   = cfertsgc(iy) * 1e-04
!            cfertsoy(iy)   = cfertsoy(iy)   * 1e-04

! temporarily assign wheat fertilizer usage to be 3 times that of
! soybeans
!            cfertwheat(iy) = cfertsoy(iy) * 3.0
!         enddo
!20       continue 
!         close(15)
!      endif

      do 100 i = lbeg, lend

         if(jday.eq.1) then
            cropout(i,1,35) = td(i)
            cropout(i,2,35) = precip(i)
            cropout(i,3,35) = xint(i,imonth)
         elseif(jday.gt.1.and.jday.lt.365) then
            cropout(i,1,35) = cropout(i,1,35) + td(i)
            cropout(i,2,35) = cropout(i,2,35) + precip(i)
            cropout(i,3,35) = cropout(i,3,35) + xint(i,imonth)
         elseif(jday.eq.365) then
            cropout(i,1,35) = (cropout(i,1,35)+td(i) ) /365.
            cropout(i,2,35) = (cropout(i,2,35)+precip(i)) / 365.
            cropout(i,3,35) = (cropout(i,3,35)+xint(i,imonth)) / 365.
         endif

         if(jday.eq.182) then
            cropout(i,4,35) = td(i)
            cropout(i,5,35) = precip(i)
            cropout(i,6,35) = xint(i,imonth)
         elseif(jday.ne.181.and.jday.ne.182.and.jday.ne.366) then
            cropout(i,4,35) = cropout(i,4,35) + td(i)
            cropout(i,5,35) = cropout(i,5,35) + precip(i)
            cropout(i,6,35) = cropout(i,6,35) + xint(i,imonth)
         elseif(jday.eq.181) then
            cropout(i,4,35) = (cropout(i,4,35)+td(i) ) /365.
            cropout(i,5,35) = (cropout(i,5,35)+precip(i)) / 365.
            cropout(i,6,35) = (cropout(i,6,35)+xint(i,imonth)) / 365.
         endif

         do 50 j = scpft, ecpft       ! crop plant functional types only
! in order to only allow a crop to be planted once each year
! initialize cropplant = 0., but hold it = 1 through the end of
! the year
!
! initialize other variables that are calculated for crops
! on an annual basis in cropresidue subroutine

            if ((j .eq. 13 .or. j .eq. 14 .or. j .eq. 16) .or. (j .eq. 15 .and. iwheattype .eq. 1)) then
!
               if (iday .eq. pdmin(j) .and. imonth .eq. pmmin(j) .and. croplive(i,j) .ne. 1.0 .and. &
                  exist(i,j) .eq. 1 .and. ncyears(i) .ge. 1) pstart(i,j)=cdays(i)

               if(iyear .gt. iyear0 .and. j .lt. 16 .and. pstart(i,j) .ne. 999 .and. &
                 (pstart(i,j)+mxmat(j)) .ge. 365) then
                     print*,j,pstart(i,j),mxmat(j),(pstart(i,j)+mxmat(j))
                     print *, 'WARNING: annual crop should not exceed end of planting year without harvest'
                     print *, 'pmmin/pdmin + mxmat should be < 365 days or before than  pcm/pcd'
                     print *, ' params.crp'
                endif

               if (iday.eq.pcd(j).and.imonth.eq.pcm(j)) then
                  cropout(i,j,1)=pdate(i,j)
                  pdate(i,j) = 0.0

                  if (croplive(i,j) .eq. 0) then
                     cropplant(i,j) = 0.0
                     if(j.eq.16) cropy(i) = 0
                  elseif (croplive(i,j) .eq. 1 .and. j .eq. 15 .and. iwheattype .eq. 2) then
                     cropplant(i,j) = 0.0
                  endif
                  harvdate(i,j)  = 999
                  cropout(i,j,2)=hdate(i,j)
                  hdate(i,j)     = 0.0
                  cropout(i,j,3)=harvidx(i,j)
                  harvidx(i,j)   = 0.0
                  cropout(i,j,4)=croplaimx(i,j)
                  croplaimx(i,j) = 0.0
                  cropout(i,j,5)=grainn(i,j)
                  grainn(i,j)    = 0.0
                  cropout(i,j,6)=cropyld(i,j)
                  cropyld(i,j)   = 0.0
                  cropout(i,j,7)=dmyield(i,j)
                  dmyield(i,j)   = 0.0
                  cropout(i,j,8)=dmleaf(i,j)
                  dmleaf(i,j)    = 0.0
                  cropout(i,j,9)=dmstem(i,j)
                  dmstem(i,j)    = 0.0
                  cropout(i,j,10)=dmroot(i,j)
                  dmroot(i,j)    = 0.0
                  cropout(i,j,11)=dmresidue(i,j)
                  dmresidue(i,j) = 0.0
                  cropout(i,j,12)=dmcrop(i,j)
                  dmcrop(i,j)    = 0.0
                  cropout(i,j,13)=residuen(i,j)
                  residuen(i,j)  = 0.0
                  cropout(i,j,14)=nconcl(i,j)
                  nconcl(i,j)    = 0.0
                  cropout(i,j,15)=nconcs(i,j)
                  nconcs(i,j)    = 0.0
                  cropout(i,j,16)=nconcr(i,j)
                  nconcr(i,j)    = 0.0
                  cropout(i,j,17)=nconcg(i,j)
                  nconcg(i,j)    = 0.0
                  cropout(i,j,18)=cropn(i,j)
                  cropn(i,j)     = 0.0
                  cropout(i,j,19)=cropfixn(i,j)
                  cropfixn(i,j)  = 0.0
                  cropout(i,j,20)=cntops(i,j)
                  cntops(i,j)    = 40.0
                  cropout(i,j,21)=cnrootvec(i,j)
                  cnrootvec(i,j)    = cnroot
                  cropout(i,j,22)=fertinput(i,j)
                  fertinput(i,j) = 0.0
                  cropout(i,j,23)=cropout(i,j,50)
                  cropout(i,j,50)=0
                  cropout(i,j,24)=crmclim(i,j)
                  cropout(i,j,25)=crmact(i,j)
                  cropout(i,j,26)=crmplant(i,j) 
                  cropout(i,j,27)=cropout(i,j,52)
                  cropout(i,j,52)=0
                  cropout(i,j,28)=idppout(i,j)
                  idppout(i,j)=0.0
                  cropout(i,j,29)=cropy(i)
                  cropout(i,j,30)=cropout(i,j,51)
                  cropout(i,j,51)=0
               endif 
            endif

            if (exist(i,j) .eq. 1. .and. croplive(i,j) .ne. 1.0 .and. cropplant(i,j) .eq. 0.0) then 
               if ((j.eq.13 .or. j.eq.14) .or. (j .eq. 15 .and. iwheattype .eq. 1)) then 
                  if (iyear .le. iyear0 + 1) then

                     if (a10td(i)  .gt. ptemp(j)    .and. &     ! 10-day average soil temperature
                        a10tmin(i) .gt. pmintemp(j) .and. &
                        cdays(i)   .ge. pstart(i,j)) then    ! impose earliest planting date
                        croplive(i,j)  = 1.0         
                        cropplant(i,j) = 1.0         
                        idop(i,j) = jday     

                        if (j .eq. 13) soydop(1,iyear-iyear0+5) = jday
                        if (j .eq. 14) corndop(1,iyear-iyear0+5) = jday
                        if (j .eq. 15) whtdop(1,iyear-iyear0+5) = jday

                        gddmaturity(i,13) = max(950., min (gdd10(i), hybgdd(j)))
                        gddmaturity(i,14) = max(950., min (gdd8(i)  * 0.90, hybgdd(j))) 
                        gddmaturity(i,15) = max(950., min (gdd0c(i), hybgdd(j)))

                     endif 
                  elseif (iyear .ge. iyear0 + 2 .and. iyear .le. iyear0 + 5) then      
                     if (a10td(i) .gt. ptemp(j) .and.     &     ! 10-day average soil temperature
                        a10tmin(i) .gt. pmintemp(j) .and. &
                        cdays(i) .ge. pstart(i,j)) then 
                        croplive(i,j) = 1.0
                        cropplant(i,j) = 1.0         
                        idop(i,j) = jday     
                        if (j .eq. 13) soydop(1,iyear-iyear0+5) = jday
                        if (j .eq. 14) corndop(1,iyear-iyear0+5) = jday
                        if (j .eq. 15) whtdop(1,iyear-iyear0+5) = jday
                        gddmaturity(i,13) = max(950., min (gddsoy(i,iyear-iyear0+5-1) *0.80, hybgdd(j)))   
                        gddmaturity(i,14) = max(950., min (gddcorn(i,iyear-iyear0+5-1) *0.90, hybgdd(j)))  
                        gddmaturity(i,15) = max(950., min (gddcorn(i,iyear-iyear0+5-1) *1.2,  hybgdd(j)))
                     endif
                  else 
                  sumdp(i,j) = 0
                  sumhy(i,j) = 0
! insert here - do iy from iyear-5 to iyear-1 --> previous 5 year mean of hybrids planted
!FIXME: gabriel abrahao: This is already implemented, but causes problems in tropical regions and is overriden right below it. See below for more information.
! keep planting date flexible for that year's weather conditions
                     yc = 5.0
                     do iy = iyear-5,iyear-1     ! hybrid based on previous 5 year average - farm management 
                        if (j .eq. 13) then 
                           sumhy(i,j) = sumhy(i,j) + gddsoy(i,iy-iyear0+5+1) * 0.8
                           sumdp(i,j) = sumdp(i,j) + (soydop(1,iy-iyear0+5+1))
                        elseif (j .eq. 14) then 
                           sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy-iyear0+5+1) * 0.9
                           sumdp(i,j) = sumdp(i,j) + (corndop(1,iy-iyear0+5+1))
                        elseif (j .eq. 15) then
                           sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy-iyear0+5+1) * 1.2
                           sumdp(i,j) = sumdp(i,j) + (whtdop(1,iy-iyear0+5+1))
                        endif
                     enddo
                     avehybrid(i,j) = sumhy(i,j) / yc
                     iavepdate(i,j) = int(sumdp(i,j)/yc)
 
                     if (a10td(i) .gt. ptemp(j) .and.     &  ! 10-day average soil temperature
                        a10tmin(i) .gt. pmintemp(j) .and. &
                        cdays(i) .ge. pstart(i,j)) then      ! impose earliest planting date 
                        
                           croplive(i,j) = 1.0
                           cropplant(i,j) = 1.0
                           idop(i,j)  = jday
                           gddmaturity(i,j) = max(950.,min(avehybrid(i,j), hybgdd(j)))

                     endif
                  endif 
!FIXME: gabriel abrahao: this basically overrides the last lines of code that change the hybrid's GDD in an adaptation strategy based on the distance between freeze events. The reason is that this obviously doesnt make sense in tropical regions where there's no freezing, and the last lines of code will basically set the hybrid to a 950 GDD one. It would be interesting to implement a different strategy here for when there is no freeze, though, probably based on the length of the rainy season on that pixel.
!override gddmaturity
gddmaturity(i,j) = hybgdd(j)

  
               endif

		if (cropy(i).eq.0) then
          if  (j.eq.16)  then 
             if (iyear .le. iyear0+1) then         
               if (a10td(i)         .gt. ptemp(j)     .and.    &     ! 10-day average soil temperature
                  a10tmin(i)       .gt. pmintemp(j)  .and. &
                  cdays(i).ge.pstart(i,j).and.cdays(i).le.(pstart(i,j)+180) )then     ! impose earliest planting date 
                  croplive(i,j)   = 1.0         
                  cropplant(i,j)  = 1.0         
                  idop(i,j)       = jday    
		  cropy(i)=1
                  sgcdop(1,iyear-iyear0+5)   = jday
                  if(gdd12(i) .eq. 0) gdd12(i)=0.000000000000001
                  if(gddsgcp(i,1) .eq. 0) gddsgcp(i,1)=0.000000000000001
                  if(gddsgcp(i,2) .eq. 0) gddsgcp(i,2)=0.000000000000001
                  gddmaturity(i,16) = max(gdd12(i)*0.8, min (gdd12(i) , hybgdd(j)))
                  gddmaturity(i,j) = gddmaturity(i,j)*(gddsgcp(i,1)/gddsgcp(i,2))

               endif 
             elseif (iyear .ge. iyear0+2 .and. iyear .le. iyear0+5) then       ! after initial spinup for crop average planting dates
               if (a10td(i)         .gt. ptemp(j)     .and.      &        ! 10-day average soil temperature
                  a10tmin(i)       .gt. pmintemp(j)  .and. &
      cdays(i).ge.pstart(i,j).and.cdays(i).le.(pstart(i,j)+180) )then     ! impose earliest planting date 
                     croplive(i,j)   = 1.0        
                     cropplant(i,j)  = 1.0         
                     idop(i,j)         = jday  
		     cropy(i)=1
                     sgcdop(1,iyear-iyear0+5)  = jday
                     if(gdd12(i) .eq. 0) gdd12(i)=0.000000000000001
                     if(gddsgcp(i,1) .eq. 0) gddsgcp(i,1)=0.000000000000001
                     if(gddsgcp(i,2) .eq. 0) gddsgcp(i,2)=0.000000000000001
                           gddmaturity(i,16) = max(gdd12(i)*0.8, min (gddsgc(i,iyear-iyear0+5-1), hybgdd(j))) ! assign hybrid based on last year
                           gddmaturity(i,j) = gddmaturity(i,j)*(gddsgcp(i,1)/gddsgcp(i,2))
 
                endif
             else 
                        sumdp(i,j) = 0
                        sumhy(i,j) = 0

!
! insert here - do iy from iyear-5 to  iyear-1 --> previous 5 year mean of hybrids planted
! keep planting date flexible for that year's weather conditions 
!
                        yc = 5.0
                        do iy = iyear-5,iyear-1     ! hybrid based on previous 5 year average - farm management 
                           if (j .eq. 16) then 
                              sumhy(i,j) = sumhy(i,j) + gddsgc(i,iy-iyear0+5+1)
                              sumdp(i,j) = sumdp(i,j) + (sgcdop(1,iy-iyear0+5+1))
                           endif
                        enddo
                        avehybrid(i,j) = sumhy(i,j) / yc
                        iavepdate(i,j) = int(sumdp(i,j)/yc)
                        if (a10td(i) .gt. ptemp(j) .and. a10tmin(i) .gt. pmintemp(j) .and.   &       ! 10-day average soil temperature
                           cdays(i) .ge. pstart(i,j) .and. cdays(i) .le. (pstart(i,j)+180)) then     ! impose earliest planting date 
                           croplive(i,j) = 1.0
                           cropplant(i,j) = 1.0
                           cropy(i) = 1
                           idop(i,j) = jday
                           if(gdd12(i) .eq. 0) gdd12(i)=0.000000000000001
                           if(gddsgcp(i,1) .eq. 0) gddsgcp(i,1)=0.000000000000001
                           if(gddsgcp(i,2) .eq. 0) gddsgcp(i,2)=0.000000000000001
                           gddmaturity(i,j) = max(gdd12(i)*0.8,min(avehybrid(i,j), hybgdd(j)))
                           gddmaturity(i,j) = gddmaturity(i,j)*(gddsgcp(i,1)/gddsgcp(i,2))

                        endif
                     endif
                  endif
               endif
!
               if (j .eq. 15 .and. iwheattype .eq. 2) then   ! plant winter wheat
!
! add check to only plant winter wheat after other crops (soybean, maize)
! have been harvested 
!
! *** remember order of planting is crucial - in terms of which crops you want
! to be grown in what order ***    
!
! in this case, corn or soybeans are assumed to be planted before
! wheat would be in any particular year that both pfts are allowed
! to grow in the same grid cell (e.g., double-cropping)  
!

                  if (iyear .eq. iyear0) then
                     if (a5tmin(i) .le. pmintemp(j) .and. imonth .ge. pmmin(j) .and. &
                        iday .ge. pdmin(j) .and. (harvdate(i,13) .ne. 999 .or.       &
                        harvdate(i,14) .ne. 999 .or. irotation .eq. 0) .and.         &
                        gdd0c(i) .ge. gddmin(j)) then
!
                        croplive(i,j) = 1.0     
                        cropplant(i,j) = 1.0     
                        idop(i,j) = jday
                        whtdop(1,iyear-iyear0+5) = jday
                        gddmaturity(i,15) = max(950., min (gdd0c(i) * 0.90, hybgdd(j)))

!
! plant winter wheat at latest possible date 
! and after all other crops were harvested for that year
!
                     elseif (imonth .ge. pmmax(j) .and. iday .ge. pdmax(j) .and.       &
                            (harvdate(i,13) .ne. 999 .or. harvdate(i,14) .ne. 999 .or. &
                            irotation .eq. 0) .and. gdd0c(i) .ge. gddmin(j)) then
                        croplive(i,j) = 1.0     
                        cropplant(i,j) = 1.0     
                        idop(i,j) = jday
!                       gddmaturity(i,15) = hybgdd(j)
                        gddmaturity(i,15) = max(950., min (gdd0c(i), hybgdd(j)))
                     endif
!
                  elseif (iyear .gt. iyear0 .and. iyear .lt. iyear0+5) then       ! after initial spinup for crop average planting dates
!
                     if (a5tmin(i) .le. pmintemp(j) .and. imonth .ge. pmmin(j) .and. &
                        iday .ge. pdmin(j) .and. (harvdate(i,13) .ne. 999 .or.       &
                        harvdate(i,14) .ne. 999 .or. irotation  .eq. 0) .and.        &
                        gdd0c(i) .ge. gddmin(j)) then
!
                        croplive(i,j) = 1.0     
                        cropplant(i,j) = 1.0     
                        idop(i,j) = jday
                        whtdop(1,iyear-iyear0+5)  = jday
                        gddmaturity(i,15) = max(950., min (gddcorn(i,iyear-iyear0+5-1) * 1.20, hybgdd(j)))  ! assign hybrid based on last year

!
! plant winter wheat at latest possible date 
! and after all other crops were harvested for that year
!
                     elseif (imonth .ge. pmmax(j) .and. iday .ge. pdmax(j) .and.        &
                            (harvdate(i,13) .ne. 999 .or. harvdate(i,14) .ne. 999 .or.  &
                            irotation .eq. 0) .and. gdd0c(i).ge. gddmin(j)) then
!
                        croplive(i,j)     = 1.0     
                        cropplant(i,j)    = 1.0     
                        idop(i,j)         = jday
                        gddmaturity(i,15) = max(950., min (gddcorn(i,iyear-iyear0+5-1) * 1.20, hybgdd(j)))  ! assign hybrid based on last year 

                     endif
!
                  else
                     sumdp(i,j) = 0
                     sumhy(i,j) = 0
!
! insert here - do iy from iyear-5 to  iyear-1 --> previous 5 year mean of hybrids planted
! keep planting date flexible for that year's weather conditions 
!
                     yc = 5.0
!                    yc = 11.0
!                    do iy = 1949, 1959          ! 11 year averaging spinup for crops
                     do iy = iyear-5,iyear-1     ! hybrid based on previous 5 year average - farm management 
                        sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy-iyear0+5+1) * 1.2
                        sumdp(i,j) = sumdp(i,j) + (whtdop(1,iy-iyear0+5+1))
                     enddo
!
                     avehybrid(i,j) = sumhy(i,j) / yc
                     iavepdate(i,j) = int(sumdp(i,j)/yc)

                     if (a5tmin(i) .le. pmintemp(j) .and. imonth .ge. pmmin(j) .and. &
                        iday .ge. pdmin(j) .and. (harvdate(i,13) .ne. 999 .or.       &
                        harvdate(i,14) .ne. 999 .or. irotation .eq. 0) .and.         &
                        gdd0c(i) .ge. gddmin(j)) then
!
                        croplive(i,j) = 1.0     
                        cropplant(i,j) = 1.0     
                        idop(i,j) = jday
                        gddmaturity(i,15) = max(950., min (avehybrid(i,j), hybgdd(j)))  ! assign hybrid based on last year 

!
! plant winter wheat at latest possible date 
! and after all other crops were harvested for that year
!
                     elseif (imonth .ge. pmmax(j) .and. iday .ge. pdmax(j) .and.       &
                            (harvdate(i,13) .ne. 999 .or. harvdate(i,14) .ne. 999 .or. &
                            irotation .eq. 0) .and. gdd0c(i).ge. gddmin(j)) then  
!
                        croplive(i,j)     = 1.0     
                        cropplant(i,j)    = 1.0     
                        idop(i,j)         = jday
!                       idop(i,j)         = iavepdate(i,j) 
                        gddmaturity(i,15) = max(950., min (avehybrid(i,j), hybgdd(j)))  ! assign hybrid based on last year 

                     endif
                  endif  ! (year check)
               endif
            endif  ! crop existence
! add fertilizer nitrogen input for each crop planted (kg m-2)
! on the planting date
! either input here, or read from a gridded dataset
! is treated a single, broadcast pulse at planting to the top soil layer
! this fertilizer is assumed to be ammonium nitrate, in a 50/50 ratio
! of NH4/NO3
!
! define amount of fertilizer added when crop is planted
! use historical changes for years between 1945-1996 (Alexander et al.,) 
! only add fertilizer on day of planting - use croplive funtion to
! make sure that in successive years, idop from previous year doesn't
! get applied here.
!
! also, don't cycle through all crop types - only the one that
! is planted this year...otherwise it will zero out the fertilizer
! values for pfts 13, 14 if going through all pfts through 15 
!
            if (jday .eq. idop(i,j) .and. croplive(i,j) .eq. 1.0 .and. exist(i,j) .eq. 1.) then
!
! constant fertilizer application each year


               fertnitro(i,13) = 0.0025      ! soybeans - kg_n m-2 y-1
               fertnitro(i,14) = 0.0180      ! maize    - kg_n m-2 y-1
               fertnitro(i,15) = 0.0080      ! wheat    - kg_n m-2 y-1
               fertnitro(i,15) = 0.0160      ! wheat    - kg_n m-2 y-1
               fertinput(i,j) = fertnitro(i,j) * 1.e+04
            elseif (exist(i,j) .eq. 1.0) then 
               fertnitro(i,j) = 0.0    
            endif                            ! planting date
            if(idpp(i,j).eq.1.and.j.eq.16) then
               fertnitro(i,16)= 0.025        !kg_n m-2 y-1 or 100 Kg/ha
            elseif (idpp(i,j).ne.1.and.j.eq.16) then
               fertnitro(i,16)= 0.0
            endif
!

50       continue
!
! end of loop
!
100   continue
!
! return to main program
! 
return
!      
end subroutine planting
