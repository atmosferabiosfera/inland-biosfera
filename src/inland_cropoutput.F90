#include "inland_config.h"
#include "inland_compar.h"
!----------------------------------------------------------------------
subroutine cropoutput(jday) 
!----------------------------------------------------------------------
!
! output crop variables for a single grid cell - diagnostic output 
!
! common blocks
!
      use inland_parameters
      use inland_control, only:iyear,imonth,iday
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
!
! local variables
!
      integer jday,i,j,k,l,n
!
      real*8 tnsoi,nuptake,hi,  &
             grain,fyld,yld
!
! define soil layer for which you want output for
! calculate total immobile/mobile N pool totals to
! specified soil depth (isoilay input)
!
      do 100 i = lbeg, lend
         tsnimm(i) = 0.
         tsnmob(i) = 0.
         tnsoi     = 0.
!
! total soil inorganic nitrogen pools available to plant for uptake
!
         do 330 k = 1, isoilay
            tsnimm(i) = tsnimm(i) + smsoil(i,k)  ! immobile pool
            tsnmob(i) = tsnmob(i) + smsoln(i,k)  ! mobile pool (e.g., nitrate) 
            tnsoi = tnsoi + smsoil(i,k) + smsoln(i,k)  
330      continue
!
! output data for model studies
! cropn is in kg/ha, totnvegn is not 
!
         nuptake = 0.0
         do 335 j = scpft, ecpft
            nuptake = nuptake + cropn(i,j)
335      continue
         nuptake = nuptake + totnvegn(i) * 1.e+04
!
! write out variables 
!
         do 380 n = scpft, ecpft
            if (exist(i,n) .eq. 1.0 .and. cropplant(i,n) .eq. 1. &
               .and. i .eq. 1 .and. harvdate(i,n) .ne. 999.) then 
!
               open(16, file = 'output/crop.output.dat',status = 'unknown')
               write(16,340) iyear,                &
                             cropyld(i,n),         &
!                            iyear,                &
!    			     n,                    &
!    			     tnsoi*1.e04,          &
                             fertinput(i,n),       &
!    			     ydeposn(i)*1.e04,     &
                             aynmintot(i)*1.e04,   &
!    			     ayimmtot(i)*1.e04,    &
                             nuptake,              &
!    			     yno3leach(i),         &
!    			     concn(i),             &
!    			     snbalance(i),         &
!    			     cropyld(i,n),         &
!    			     dmyield(i,n),         &
                             harvidx(i,n),         &
!    			     dmleaf(i,n),          &
!    			     dmstem(i,n),          &
!    			     dmroot(i,n),          &
!    			     dmyield(i,n),         &
!    			     dmcrop(i,n),          &
!    			     cntops(i,n),          &
!    			     cnroot(i,n),          &
                             nconcl(i,n),          &
                             nconcs(i,n),          &
                             nconcr(i,n),          &
                             nconcg(i,n),          &
                             grainn(i,n),          &
!    			     drntot(i),            &
!    			     ayprcp(i),            &
                             croplaimx(i,n),       &
                             cropfixn(i,n)        
!    			     totcsoi(i),           &
!                             idop(i,n),            &
!                             harvdate(i,n)         
            endif
!
340            format(i6,5f8.2,4f8.4,3f8.2,i5,i5)
!340           format(i6,f8.2)
!
380      continue  ! loop through all pfts
100   continue
!
      return
!
end subroutine cropoutput
