#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine cropresidue(jday)
! ---------------------------------------------------------------------
!
! routine that calculates the amount of residue that is input back to
! the soil at the end of growing season from crop at harvest
!
! also keeps track of effects of residue and plant nitrogen uptake on
! the nitrogen budget
!
!
! common blocks
!
      use inland_parameters
      use inland_control, only:iyear,iyear0, imonth
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
!
!
      real*8 dumg,             &
!	     maxhi(npft),      &   ! Maximum harvest index allowed for crops
!   	     fyield(npft),         ! Adjustment factor for fraction of cbiog that 
                                   ! is just grain - measured in the field 
!   	     cgrain(npft),     &   ! carbon fraction in grain
!   	     convfact(npft),   &   ! conversion factor to go to bu/acre 
             rdm,              &
             fdml,             &
             fdms,             &
             sumhy(npoi,npft), &
             fnresidue, yc
!
      integer jday, i, j, iy
!
!
! begin global grid
!
      do 100 i = lbeg, lend
!
         if (isimagro .gt. 0) then
!
! zero out litter fall rates
!  
         falll(i) = 0.0
         fallw(i) = 0.0 
         fallr(i) = 0.0   
!
            if(croplive(i,16).eq. 1) then
              fallrsgc(i,2) = max(nint(fallrsgc(i,2) - fallrsgc(i,1)*(1./90.)), 0)
                 if(fallrsgc(i,2) .gt. 0) fallr(i) = fallr(i) + fallrsgc(i,1) * (1./90.)
                   fallr(i) = fallr(i) + fallrsgc(i,3)
            endif
!
            do 90 j = scpft, ecpft
!
! calculate CRM values from Pioneer regression relationships 
! Jan 28 03 CJK
!
               crmclim(i,j)    = max(73., min((gddmaturity(i,j)+53.683)/13.882,135.))
               crmact(i,14)    = max(73., min((gddfzcorn(i)+53.683)/13.882,135.))
               crmact(i,16)    = max(73., min((gddfzsgc(i)+53.683)/13.882,135.))

!
! gddplant gets reinitialized to 0.0 at maturity date so save value here
!
               if (gddplant(i,j) .gt. 0.0 .and. croplive(i,j) .eq. 1.) then
                  crmplant(i,j)   = max(73., min((gddplant(i,j)+53.683)/13.882,135.))
               endif
!
! only write out values at harvest date, and re-initialize crop variables
! at this time - this allows for the same crop (e.g., wheat) to be grown
! across two consecutive calendar years   
!	
               if (exist(i,j) .eq. 1.0 .and. harvdate(i,j) .eq. jday) then
                  pdate(i,j) = idop(i,j)
                  idppout(i,j) = idpp(i,j)
                  hdate(i,j) = harvdate(i,j)
                  ayabprod(i,j) = max(cbiol(i,j), ayabprod(i,j))
!
! added so division by zero doesn't take place in those places where
! production was zero 
!
                  if(j.ne.16) then
                     dumg       = cbiog(i,j) 
                     cbiog(i,j) = cbiog(i,j) * fyield(j) 

                     cbios(i,j) = max(0.0, cbios(i,j) + (dumg - cbiog(i,j)))
!
! impose upper limit on harvest index as per JMN suggestion 
! if harvest index is > allowable, put excess to stem which will add to more
! litterfall and adjust the grain accordingly
! might have to revisit logic with respect to what actually get harvested, versus
! what is left in the field after harvest as litter input to the soil 
! 
                     harvidx(i,j)   =  cbiog(i,j) / ayabprod(i,j)
                     croplaimx(i,j) = plaimx(i,j)   
!
                     if (harvidx(i,j) .gt. maxhi(j)) then
                        harvidx(i,j) = maxhi(j) 
                        dumg         = cbiog(i,j) 
                        cbiog(i,j)   = maxhi(j) * ayabprod(i,j)
! 
! add excess to stem of plant
!
                        cbios(i,j)   = cbios(i,j) + (dumg - cbiog(i,j))
!
                     endif

                  endif 
!
                  grainn(i,j)   = (cbiog(i,j) / cfrac(j)) * fngrain(i,j) * 1e+04
!
                  if(j.eq.16)then
                     cropyld(i,j)=(cbiog(i,j)+cbios(i,j))*fyield(j)*(1./0.3)*10.0*(1./cgrain(j))*1.07 
                  else           
                     cropyld(i,j) = cbiog(i,j) * convfact(j) / cgrain(j)
                  endif
!
                  dmyield(i,j) = cbiog(i,j)* 10.0 /cgrain(j)  
                  dmstem(i,j)  = cbios(i,j)* 10.0 /cgrain(j)   
                  dmleaf(i,j)  = cbiol(i,j)* 10.0 /cgrain(j)   
                  dmroot(i,j)  = cbior(i,j)* 10.0 /cgrain(j)   
                  dmcrop(i,j)  = dmyield(i,j) + dmstem(i,j) + dmleaf(i,j) + dmroot(i,j)
!
                  if(j .eq. 16)then
                     dmyield(i,j) = dmyield(i,j) * fyield(j) 
                     dmstem(i,j)  = dmstem(i,j)  * fyield(j)
                  endif
!
                  if(j .eq. 16) then             
                     dmresidue(i,j) = dmleaf(i,j)+( (cbios(i,j)+cbiog(i,j))*(1-fyield(j))*10.0/cgrain(j) ) 
                  else
                     dmresidue(i,j) = dmleaf(i,j)  + dmstem(i,j)
                  endif
!
                  rdm  = dmresidue(i,j) + (aylprod(i,j)*10.0 /cgrain(j)) - dmleaf(i,j)
                  fdml = (aylprod(i,j)*10.0 /cgrain(j)) / rdm
                  fdms = (dmresidue(i,j)-dmleaf(i,j)) / rdm
                  fnresidue  = fnleaf(i,j) * fdml + fnstem(i,j) * fdms  

!
! calculate amount of N in aboveground residue (kg/ha) 
!
                  residuen(i,j)  = (fnleaf(i,j) * dmleaf(i,j) + fnstem(i,j) * &
                                   (dmresidue(i,j)-dmleaf(i,j))  ) * 1e+04
!
! assign leaf, stem, root, and grain nitrogen concentrations
! to new variables (in percent)
!
                  nconcl(i,j)    = fnleaf(i,j)   * 100.0
                  nconcs(i,j)    = fnstem(i,j)   * 100.0
                  nconcr(i,j)    = fnroot(i,j)   * 100.0
                  nconcg(i,j)    = fngrain(i,j)  * 100.0
!   
! assign total nitrogen plant uptake to new variable for
! purposes of outputting data at end of year (kg/ha)
!
                  cropn(i,j)     = totnuptake(i,j) * 1.e+04  
!
! assign total nitrogen fixation for crop to new variable for
! purposes of outputting this data at end of calendar year
! units of kg/ha
!
                  cropfixn(i,j)  = totnfix(i,j) * 1.e+04
!
! carbon nitrogen ratio of plant residue goes to biogeochem.f
! and fine roots
! 

                  if (fnresidue .gt. 0.0) then
                     cntops(i,j) = min(cfrac(j) / fnresidue, 200.0)
                  else
                     cntops(i,j) = 60.0
                  endif
!
                  if (fnroot(i,j) .gt. 0.0) then
                     cnrootvec(i,j) = min(cfrac(j) / fnroot(i,j), 200.0) 
                  else
                     cnrootvec(i,j) = 80.0
                  endif
!  
! assume that stem and leaf are both included in leaf litterfall value
! these annual total values (falll, fallr, fallw) cannot be changed 
! on a daily basis because biogeochem.f uses the annual total, split
! between each day of the year equally 
! carbon returned as residue
!
                  if (j.eq.16) then
                  
                     if(firecane .eq. 1) then
                        falll(i) = falll(i) + (cbios(i,j)+cbiog(i,j))*(1-fyield(j))/2.  
                     else
                        falll(i) = falll(i) + aylprod(i,j) + (cbios(i,j)+cbiog(i,j))*(1-fyield(j)) 
                     endif
                  else
                     falll(i) = falll(i) + aylprod(i,j) + cbios(i,j) 
                  endif
!
                  if (j .eq. 16 .and. cropy(i) .le. nratoon) then
                     fallr(i) = fallr(i) + cbior(i,j) * 0.30
                     fallrsgc(i,1) =  cbior(i,j) * 0.70
                     fallrsgc(i,2) =  cbior(i,j) * 0.70

                  elseif (j .eq. 16) then
                     fallr(i) = fallr(i) + cbior(i,j) 
                     fallrsgc(i,1) =  cbior(i,j) * 0.
                     fallrsgc(i,2) =  cbior(i,j) * 0.
                  else
                     fallr(i) = fallr(i) +  ayrprod(i,j) 
                  endif
!
                  fallw(i)       = fallw(i) +  cbiow(i,j)
                  plai(i,j)      = 0.01
                  thrlai(i,j)    = 0.0
                  peaklai(i,j)   = 0.0 - rwood * cbiow(i,1) * sapfrac(i) * funca - &
                                         rroot * cbior(i,1) * funcb
                  ccdays(i,j)    = 0.0
                  cbiol(i,j)     = 0.0
                  cbior(i,j)     = 0.0        
                  cbios(i,j)     = 0.0
                  cbiog(i,j)     = 0.0
                  cbiow(i,j)     = 0.0
                  hui(i,j)       = 0.0
                  aybprod(i,j)   = 0.0
                  ayrprod(i,j)   = 0.0
                  ayabprod(i,j)  = 0.0
                  aylprod(i,j)   = 0.0
                  leafout(i,j)   = 0.0
                  htmx(i,1)      = 0.0
                  cumlvs(i,j)    = 0.0
                  plaimx(i,j)    = 0.0
                  dpgf(i,j)      = 0.0
                  biomass(i,j)   = 0.0
                  totnuptake(i,j)= 0.0
                  tnplant(i,j)   = 0.0
                  totnfix(i,j)   = 0.0
                  idpp(i,j)      = 0.0
                  cropout(i,j,50)= gddmaturity(i,j)
                  cropout(i,j,51)= gddplant(i,j)
                  cropout(i,j,52)= grainday(i,j)
                  gddplant(i,j)  = 0.0
                  gddtsoi(i,j)   = 0.0
                  grainday(i,j)  = 9999.
                  sai(i,1)       = 0.0
                  fu(i)          = 0.0
                  lai(i,1)       = 0.0
                  zbot(i,1)      = 0.0
                  ztop(i,1)      = 0.0
                  totbiol(i)     = 0.0
                  totlail(i)     = 0.0  
                  vf(i)          = 0.0  ! vernalization factor for winter wheat
                  arepr(j)       = 0.0
                  idop(i,j)      = 999

                  if(j .eq. 16 .and. cropy(i) .eq. 1) then
                  else
                     if(iyear .le. iyear0+5) then
                        gddsgcp(i,2)=max(0.0000001, (gddsgcp(i,2)+cropout(i,j,50))/2.)
                     else
                        gddsgcp(i,2)=max(0.0000001, (gddsgcp(i,2)*(nratoon-1)+cropout(i,j,50))/nratoon)
                     endif
                  endif

                  if(j .eq. 16 .and. cropy(i) .le. nratoon) then
                     cropy(i)      = cropy(i)+1
                     croplive(i,j) = 1.0     
                     idop(i,j)     = jday+1                    
!
                     if(ccdays(i,j) .ge. 1 .and. idppout(i,j) .le. mxmat(j)-30) then
                        croplive(i,j)     = 0.    
                        idop(i,j)         = 999.
                     endif
!
                     if(ccdays(i,j) .ge. 1 .and. cropy(i) .eq. 2) then   
                        croplive(i,j)     = 0.     
                        idop(i,j)         = 999.
                     endif
!
                     if (iyear .le. iyear0+1) then
                        gddmaturity(i,16) = max(2000., min (gdd12(i) , hybgdd(j)))
                        elseif (iyear .ge. iyear0+2 .and. iyear .le. iyear0+5) then
                        gddmaturity(i,16) = max(2000., min (gddsgc(i,iyear-iyear0+5-1), hybgdd(j))) 
                     else 
                        sumhy(i,j) = 0
                        yc = 5.0
                        do iy = iyear-5,iyear-1     
                           sumhy(i,j) = sumhy(i,j) + gddsgc(i,iy-iyear0+5+1) 
                        enddo
                        avehybrid(i,j)    = sumhy(i,j) / yc
                        gddmaturity(i,j) = max(2000.,min(avehybrid(i,j), hybgdd(j))) 
                     endif
                  else
                     cropy(i) = 0
                  endif  

               endif 

90          continue
!
         endif 
!  
100   continue       
!
! call to map out vegetation classes
! subroutine is within vegetation.f
!
      call vegmap
!
! return to main program
!
      return
!      
end subroutine cropresidue
