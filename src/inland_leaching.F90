#include "inland_config.h"
! last modified C. Kucharik 10.02.02
!--------------------------------------------------------   
subroutine leaching(irestart, iyrrestart,istep,iholdsoiln) 
!--------------------------------------------------------
! CJK November 2000
!
! algorithms were derived from the ALEXI model (Anderson et al., 1999)
! it keeps track of inorganic nitrogen pools and its location
! and movement through the soil column profile.  Each layer has a total
! quantity of solute - where all is available to the plant for uptake -
! however, a buffering constant keeps a fixed % N in soil solution at any time while the rest
! is assumed bound to soil aggregates.  This allows for a buffering action to
! be built in where as solution N is decreased, N is allowed to move from
! the total soil pool to solution. 
!
! the inputs of nitrate to the model are : fertilizer, nitrogen deposition
! nitrogen fixation in soybean, nitrogen mineralization through biogeochemistry
!
! the outputs are the leaching through the bottom of the soil profile, and
! plant n uptake through carbon assimilation
!
!
! common blocks
!
      use inland_parameters
      use inland_control,only:iday, imonth, iyear, iyear0
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
!
! local variables
!    
      integer istep,      &
              irestart,   &
              iyrrestart,   &
              iholdsoiln, &
              i,j,m,      &
              k,l,n
!
      real*8  tprod, &   ! carbon assimilation for timestep by natural vegetation (mol co2 m-2)
!              co,    &   ! initial total solute mass for the profile (kg-solute m-2) 
!              cf,    &   ! buffering factor between inorganic N in soil vs. solution 
              cassn, &
              bal1,  &   ! solute balance calculation
              bal3,  &   ! solute balance calculation
              bal,   &   ! solute balance calculation after check 
              dif        ! difference calculation resulting from balance   
!
      real*8  cnwood,    &  ! c/n ratio of wood       - same value as in biogeochem.f
!              cnfroot,   &  ! c/n ratio of fine roots - same value as in biogeochem.f
!              cnleaf,    &  ! c/n ratio of leaf	   - same value as in biogeochem.f 
!              cndepth,   &
              frootavg,  &
              ypnuptake, &
              grain,     &
              fyld,      &
              yld,       &
              hi,        &
              tnsoi,     &
              deficit,   &
              sum,       &
              rwork,     &
              rwork2,    &
              sumcrop,   &
              fnuptake,  &
!              nloss,     &
              tupsoil
!              fnitrate     ! fraction of fertilizer addition that is nitrate, and will
                           ! be kept track of for leaching and concentration in solution
                           ! fraction by weight 
!
      real*8 depth(nsoilay), &  ! cumulative depth to bottom of each soil layer
             snode(nsoilay)
!
      real*8 fin(lbeg:lend, nsoilay), &   ! nitrogen input to layer 
             drn(lbeg:lend,nsoilay+1)  ! drainage out of layer (mm/timestep)

!
! set constants
!
      rwork  = dtime * 12.e-03
!
! these defined variables are also in biogeochem.f for decomposition
! to maintain proper relationship between nitrogen mineralization and
! nitrogen uptake here due to vegetation growth requirements, make sure
! these values are consistent for natural vegetation
!
!>>>This variables are transferred for soil params and are being read in readparam.  
!
!      cnwood  = 200.0
!      cnfroot =  80.0
!      cnleaf  =  60.0
!
!      co = 0.0050      ! initial inorganic nitrogen in soil storage pools 
!      cf = 0.07        ! buffering constant for inorganic N in soil solution
!      fnitrate = 0.90  ! fraction of leached inorganic nitrogen that is nitrate - other portion is ammonium
                       ! values from kr brye
!      nloss   = 0.50   ! fraction of immobile inorganic nitrogen that is lost annually to volatilization/denitrification 

!      cndepth = 0.0
! 
! set day of year counter
!
      if (iday .eq. 1 .and. imonth .eq. 1 .and. istep .eq. 1) idoy = 0 
         if (istep .eq. 1) idoy = idoy + 1
!------------------------------------------------------------------------------
! calculate cumulative depth to bottom of each layer
! this only needs to be called once per run
! depth could be local variable, where snode will be global
!------------------------------------------------------------------------------
!
         do 110 k = 1, nsoilay 
            if (k.eq.1) then
               depth(k) = hsoi(k)
            else
               depth(k) = depth(k-1) + hsoi(k)
            endif
 110     continue
!
         do 120 k = 1, nsoilay
            if (k.eq.1) then 
!              snode(k) = (depth(k+1) - depth(k)) / 2.0 
               snode(k) = (depth(k+1) - 0.) / 2.0 
            else if (k.eq.nsoilay) then
               snode(k) = (depth(k) - depth(k-1)) / 2.0 
            else
               snode(k) = (depth(k+1) - depth(k-1)) / 2.0          
            endif
               cndepth = cndepth + snode(k) 
 120     continue
!
!------------------------------------------------------------------------------
! begin grid
!
            do 100 i = lbeg, lend  
!        
! set concentrations for first year of restart (beginning)
!
            if (irestart .eq. 1 .and. iyear .eq. iyrrestart .and. &
                istep .eq. 1 .and. iday .eq. 1 .and. imonth .eq. 1) then
               ctoti(i) = 0.0
!	      
               do 105 k = 1, nsoilay
                  if (iholdsoiln .eq. 0) then
                     smsoil(i,k) = smsoil(i,k) * 0.0
                     smsoln(i,k) = smsoln(i,k) * 0.0

                  endif
                  csoln(i,k)  = smsoln(i,k) * (1.e+09)/ &
                               (1.e+04 * (wsoi(i,k) + & 
                               wisoi(i,k)) * snode(k) * 10.)   
                  ctoti(i)    = ctoti(i) + smsoil(i,k) + smsoln(i,k)
105            continue
               ctot(i)     = ctoti(i)
            endif
!      
! set initial concentrations for first timestep of model run
!

            if (iyear.eq.iyear0 .and. istep.eq.1 .and. &
                iday.eq.1 .and. imonth.eq.1) then 
!

               ctoti(i)    = 0.0  ! total initial inorganic nitrogen in profile
               ctot(i)     = 0.0  ! total current inorganic nitrogen in soil profile (solution and soil) (kgm-2)
               drntot(i)   = 0.0  ! total drainage 
               ftot(i)     = 0.0  ! total leaching of inorganic  nitrogen
               yno3leach(i)= 0.0  ! total annual nitrate leaching
               assimn(i)   = 0.0   
               snbalance(i)= 0.0  ! soil nitrogen balance
!
               do 130 k = 1, nsoilay
!                 smsoil(i,k) = co * (snode(k) / depth(nsoilay))
                  smsoil(i,k) = co * (snode(k) / cndepth)
                  smsoln(i,k) = cf * smsoil(i,k)
!

! adjust mass balance
!
                  smsoil(i,k) = smsoil(i,k) - smsoln(i,k)
!
! check conversion factor from kg m-2 (ibis) to kg ha-1 (kris' equations) to mg liter-1
! concentration in solution
! have to account for ice fraction because wsoi can be 0.0 for much of the year 
!
                  csoln(i,k) = smsoln(i,k) * (1.e+09)/ &
                               (1.e+04 * (wsoi(i,k) + wisoi(i,k)) * snode(k) * 10.)
                  ctoti(i)    = ctoti(i) + smsoil(i,k) + smsoln(i,k)
130            continue
!      
            else if (iday.eq.1 .and. imonth.eq.1 .and. istep.eq.1) then
!
! initialize annual drainage and solute leached to zero at beginning of year
!
               drntot(i)   = 0.0     ! total drainage through profile (mm)
               ftot(i)     = 0.0     ! cumulative inorganic nitrogen leached (kg m-2)
               yno3leach(i)= 0.0     ! annual total nitrate leaching (kg no3 m-2)
               assimn(i)   = 0.0
               totnvegn(i) = 0.0     ! cumulative inorganic nitrogen uptake by natural vegetation (kg m-2 y-1)  
               taninp(i)   = 0.0     ! total annual inputs of inorganic nitrogen 
               snbalance(i)= 0.0     ! annual soil nitrogen balance
!
               ctoti(i)     = ctot(i) ! initial total inorganic N in profile at beginning of timestep 
               tsinp(i)     = 0.0     ! total daily inorganic N input
               tslay(i)     = 0.0     ! total daily inorganic N input to depth determined by soilay
               dtnleach(i)  = 0.0     ! daily total inorganic nitrogen solute (ammonium and nitrate) leached
               dnileach(i)  = 0.0     ! daily total nitrate-nitrogen leached
               tpnuptake(i) = 0.0     ! total daily plant nitrogen uptake
               ctot(i)      = 0.0     
               ddrn(i)      = 0.0     ! daily total drainage at specified depth in profile 
!
            else if (istep .eq. 1) then

!
! initialize daily sums used in mass balance calculation
!
               ctoti(i)     = ctot(i) ! initial total inorganic N in profile at beginning of timestep 
               tsinp(i)     = 0.0     ! total daily inorganic N input
               tslay(i)     = 0.0     ! total daily inorganic N input to depth determined by soilay
               dtnleach(i)  = 0.0     ! daily total inorganic nitrogen solute (ammonium and nitrate) leached
               dnileach(i)  = 0.0     ! daily total nitrate-nitrogen leached
               tpnuptake(i) = 0.0     ! total daily plant nitrogen uptake
               ctot(i)      = 0.0     
               ddrn(i)      = 0.0     ! daily total drainage at specified depth in profile 
!
! add denitrification/volatilization loss on daily timestep
!
!              do 125 k = 1, nsoilay
!                 smsoil(i,k) = smsoil(i,k) * (1. - (nloss/ndaypy))
!
!125           continue
            else  ! all other timesteps
               fin(i,1)  = 0.0     ! initialize input to top layer each timestep
               ctot(i)   = 0.0
            endif
!
!------------------------------------------------------------------------------
! calculating drainage (mm) through each node per timestep
! note that the drainage is from the above layer through the top of (k) - thus
! to get drainage out of the profile, have to assume that an additional hypothetical
! layer exists below, where bperm influences the gdrain calculated in soil.f  
!------------------------------------------------------------------------------
!
!  
!
            drn(i,1) = max(dble(0.), (fwtop(i) + fwpud(i))* dtime)  ! water infiltration (mm) into top layer
            do 140 k = 2, nsoilay+1      ! drainage out the bottom is an additional soil layer
               drn(i,k) = max(dble(0.), wflo(i,k) * dtime)  ! rate of drainage through each node -
!                                                      ! calculated in soil.f (kg m-2 s-1)
!                                                      ! drn(i,nsoilay+1) is equal to gdrain(i)

140         continue
!
! designate which layer in the profile you are tracking drainage 
! i.e. KRB measures drainage at ARL at 1.4 m with lysimeter
!
!           drntot(i) = drntot(i) + drn(i,nsoilay+1) 
            drntot(i) = drntot(i) + drn(i,isoilay) 
            ddrn(i)   = ddrn(i)   + drn(i,isoilay)
!
! determine nitrogen uptake by natural vegetation
!
!
             tprod = 0.0
            do 215 j = 1, 12     ! loop through natural vegetation types to check for carbon
                                 ! assimilation 
!
! total carbon assimilation for the timestep (tnpp has units of mol-co2 /m-2 / s)

!
              if (tnpp(i,j) .eq. -0.) tnpp(i,j) = 0.
               tprod = tprod + tnpp(i,j) * dtime
215         continue
!
! if natural vegetation has positive carbon assimilation for this timestep then calculate 
! how much nitrogen was required for growth based on carbon/nitrogen ratios dictated in biogeochem.f
!
             cassn = 0.0
            if (tprod .gt. 0.) then 
               do 220 j = 1, 12    ! loop through natural vegetation types
!
! calculate assimilated carbon for this timestep and nitrogen required
! to satisfy growth
!
                  cassn = cassn + (tnpp(i,j) * rwork * awood(i,j) / cnwood) + &
                                  (tnpp(i,j) * rwork * aleaf(i,j) / cnleaf) + &
                                  (tnpp(i,j) * rwork * aroot(i,j) / cnroot)
!
220            continue
!
               assimn(i) = assimn(i) + cassn
               tupsoil = 0.0 
!
! calculate total transpiration stream for both [u]pper and [l]ower
! canopies - calculated in canopy.f
!
               do 230 l = 1, nsoilay
                  tupsoil = tupsoil + (upsoiu(i,l) + upsoil(i,l)) * dtime 
230            continue
!
! determine that the fraction of the total nitrogen required for uptake
! is proportional to the fraction of the total transpiration stream
! contribution of that layer 
!
! anuptake is same variable used in crops.f for actual nitrogen uptake for that
! timestep
!
               do 240 l = 1, nsoilay
                  if (tupsoil .gt. 0.0) then
                     fnuptake = ((upsoiu(i,l) + upsoil(i,l)) * dtime) / tupsoil 
                  else
                     fnuptake = 0.0
                  endif
!
! add total nitrogen uptake for whole grid cell - could eventually
! be a combinatation of natural vegetation and crops
! anuptake is also being calculated in crops.f 
! 
                  anuptake(i,l) = min(dble((smsoil(i,l)+smsoln(i,l))), fnuptake * cassn)
!
! update annual nitrogen uptake by natural vegetation (forests/grasses/shrubs) 
!
                  totnvegn(i)   = totnvegn(i) + anuptake(i,l) 
240            continue
!
            endif  ! natural vegetation nitrogen uptake 
!
!------------------------------------------------------------------------------
!
            do 200 k = 1, nsoilay
!
! calculate the average root profile in the grid cell for upper and lower canopies
! combined
!
               frootavg = (froot(k,1) + froot(k,2)) / 2.0

               if (k .eq. 1 .and. istep .eq. 1) then 
!
! first timestep/top layer - each day checked
!
                  do 210 j = 1, npft
!               
! managed crop ecosystems: 
! human dependency - get at time of planting
! this could be changed in crops.f - to change timing/method of nitrogen fertilizer management
! value of fertnitro(j) to be 0 at other times of year
!
                     if (exist(i,j) .ne. 0.) then
                        if (fertnitro(i,j) .gt. 0.0 .and. croplive(i,j) .eq. 1.) then
                            fin(i,1) = fertnitro(i,j)
                        else 
                            fin(i,1) = 0.
                        endif 
                     endif
210               continue
!
! add nitrogen deposition to the amount coming into the top
! layer for that particular day - based on daily precipitation
! and calculated in biogeochem.f 
! only added on beginning timestep of each day
!
                  fin(i,1) = fin(i,1) + deposn(i)

               else      
                  fout(i,0) = 0.0
                  nout(i,0) = 0.0
                  fin(i,k)  = fout(i,k-1) 
               endif
!
! nitrogen movement - potential inputs are fertilizer, n-deposition, n-fixation,n-mineralization.  
! nitrogen mineralization is calculated in biogeochem.f as a daily rate and
! was converted to mole-N s-1
!
! have to assume where nitrogen mineralization is taking place in the profile -
! and how amount of atmopheric n-fixation (daily quantity) is being added by roots
! during this timestep (dtime)
! use a weighted average of fine root distribution profiles for lower and upper plant canopies
! fine root distribution is froot - initialized in initial.f  
!
! value for anuptake is calculated in nitrostress - which is applied to
! the vmax rate for each timestep 
!
! if crops are planted, no n-fixation used from biogeochem.f 
! that n-fixation is assumed from natural vegetation types that fix atmospheric N2
!
! n-fixation from soybeans is incorporated into soil inorganic N pools in crops.f
! 
!              if (icropsum(i) .gt. 0.) then
!                 fixsoin(i) = 0.0
!                 yfixsoin(i) = 0.0
!              endif
!
! assume natural nitrogen fixation = 0 
! cjk 11.18.01
!
               fixsoin(i)  = 0.0
               yfixsoin(i) = 0.0
               smsoil(i,k) = smsoil(i,k) + fin(i,k) + &
                             fixsoin(i) * frootavg * dtime / 86400. + &
                             (tnmin(i) * frootavg * dtime * 0.014) - &
                             anuptake(i,k)
!
! if total plant nitrogen uptake cannot be derived entirely from the immobile pool
! then remove it from the solution pool
!
               if (smsoil(i,k) .lt. 0.0) then
                  deficit     = smsoil(i,k) 
                  smsoln(i,k) = max(dble(0.), smsoln(i,k) + deficit) 
                  smsoil(i,k) = 0.0
               endif
!
! keep sum of daily plant nitrogen uptake for mass balance calculation
!
               tpnuptake(i) = tpnuptake(i) + anuptake(i,k) 
               smsoil(i,k) = smsoil(i,k) - (cf*smsoil(i,k) - smsoln(i,k))
               smsoln(i,k) = smsoln(i,k) + (cf*smsoil(i,k) - smsoln(i,k))
!
! account for ice fraction when calculating concentration in solution
!
               csoln(i,k)  = smsoln(i,k) * (1.e+09)/ &
                            (1.e+04 * (wisoi(i,k) + wsoi(i,k)) * snode(k) * 10.)
!
! drainage is designated by subscript of layer it is going INTO 
! fout is for current layer 
!
! KRBs units are kg ha-1  - for fout
! since we are kg m-2     - need to reduce fout by larger constant
!
! we are interested in tracking nitrate only for leaching - add factor in these
! equations to account for the fact that csoln can contain both ammonium and nitrate 
! the fertnitro addition assumes that both nitrate and ammonium are added together
!
!              fout(i,k)   = min(smsoln(i,k), csoln(i,k) * drn(i,k+1) / 100.0)
               fout(i,k)   = min(smsoln(i,k), csoln(i,k) * &
                             drn(i,k+1) / 1e+06)
               nout(i,k)   = min(smsoln(i,k), csoln(i,k) * fnitrate * &
                             drn(i,k+1) / 1e+06)
               fout(i,k)   = max(dble(0.), fout(i,k))
               nout(i,k)   = max(dble(0.), nout(i,k))
!
! remove leached inorganic-N from layer reassign total nitrogen to each layer available to leach/total
! redistribute inorganic-N in layer between solution and total in soil
!
               smsoln(i,k) = smsoln(i,k) - fout(i,k)
               smsoil(i,k) = smsoil(i,k) - (cf * smsoil(i,k) - smsoln(i,k))
               smsoln(i,k) = smsoln(i,k) + (cf * smsoil(i,k) - smsoln(i,k))

               if (smsoil(i,k) .lt. 0.0) then
                  smsoil(i,k) = 0.0
                  smsoln(i,k) = 0.0
               endif
               ctot(i) = ctot(i) + smsoln(i,k) + smsoil(i,k)
!
! account for ice fraction in soil when calculating solute concentration
! 
               csoln(i,k)  = smsoln(i,k) * (1.e+09)/ (1.e+04 * &
                             (wisoi(i,k)+ wsoi(i,k)) * snode(k) * 10.)
!
! assign daily nitrate concentration for 1.4 m layer and
! total amount of daily drainage through that layer for arlington comparisons 
!
               if (istep .eq. 86400. / dtime) then
                  daynconc(1,idoy) = fnitrate * adcsoln(i,isoilay)
                  daydrn(1,idoy)   = ddrn(i)
               endif
!             
200         continue  ! loop through all soil layers
!
! total inputs into the soil for this daily timestep
!
            tsinp(i) = tsinp(i) +  fin(i,1) + fixsoin(i) * dtime / 86400. + &
                       tnmin(i) * dtime * 0.014
!
! calculate the total daily inputs to top number of layers - determined by the
! soilay constant - to help in mass balance calculation to that depth
!
            do 280 k = 1, isoilay
               frootavg    = (froot(k,1) + froot(k,2)) / 2.0
               tslay(i)    = tslay(i) + tnmin(i) * frootavg * dtime * 0.014 + &
                            fixsoin(i) * frootavg *  dtime / 86400.
280         continue
!
            tslay(i) = tslay(i) + fin(i,1)
!
! calculate total amount of total nitrogen leached out of entire profile
! for daily timestep
!
            dtnleach(i) = dtnleach(i) + fout(i,nsoilay)
!
! convert to a rate (y-1) at end of day for daily output for Simon
! base on rate kg nitrate per hectare - for layer we are interested in designated
! as input to baseflow
!
            dnileach(i) = dnileach(i) + nout(i,isoilay)
!             
! update annual total nitrate-nitrogen leaching - covert to kg/ha
! trying to compare to KRBs measurements at 1.4 m in profile
! or input to baseflow at soilay = 1.5 m for regional modeling 
!
            ftot(i)      = ftot(i)      + fout(i,isoilay)  * 1e+04
            yno3leach(i) = yno3leach(i) + nout(i,isoilay)  * 1e+04
!           ftot(i)      = ftot(i)      + nout(i,nsoilay) * 1e+04
!           ftot(i)      = ftot(i)      + nout(i,isoilay)  * 1e+04
!
! end of year calculation for flow-weighted mean nitrate concentration 
!
            if (istep .eq. 86400. / dtime .and. &
               imonth .eq. 12 .and. iday .eq. 31) then 
!
! put in check for division by zero for annual drainage
! 

               if (drntot(i) .le. 0) then
                   concn(i) = 0.0
               else
                  sum = 0
                  do 250 l = 1, idoy
                     sum = sum + (daydrn(1,l)/drntot(i)) * daynconc(1,l)
250               continue
                  concn(i) = sum
               endif
            endif
!
!------------------------------------------------------------------------------
! calculate mass balance approach at each day to make sure solute is being
! conserved 
!------------------------------------------------------------------------------
!
! calculate each timestep addition of nitrogen fixation from nitrostress routine
! in crops.f 
!
            do 275 j = 1, npft
!
! only add fixed nitrogen to taninp for the top soil layers according to soilay
!
               do 290 k = 1, isoilay
                  taninp(i) = taninp(i) + fixn(i,j) * froot(k,1)
290            continue
!
! add fixation - which is a total for each timestep - to the total soil inputs
! 
               tsinp(i)  = tsinp(i)  + fixn(i,j)
!
275         continue
!
            if (istep .eq. (86400. / dtime)) then
              taninp(i) = taninp(i) + tslay(i) 
               bal1 = ctoti(i) + tsinp(i) - &
                      dtnleach(i) - ctot(i) - &
                      tpnuptake(i)
!
               if (bal1 .ne. 0.0) then      ! excess inputs vs. outputs 
                  dif = bal1
                  ctot(i) = 0.0
                  do 320 k = 1, nsoilay
!                    smsoil(i,k) = smsoil(i,k) + (dif*snode(k) / depth(nsoilay)) 
                     smsoil(i,k) = smsoil(i,k) + (dif*snode(k) / cndepth) 
                     smsoil(i,k) = smsoil(i,k) - (cf * smsoil(i,k) - smsoln(i,k))
                     smsoln(i,k) = smsoln(i,k) + (cf * smsoil(i,k) - smsoln(i,k))
                     ctot(i)     = ctot(i) + smsoil(i,k) + smsoln(i,k)
320               continue
               endif
!
! recalculate the solute balance
!
!              bal = ctoti(i) + tsinp(i) - dtnleach(i) - &
!                    ctot(i)  - tpnuptake(i) 
!
! convert to a rate (y-1) at end of day for daily output for Simon
! base on rate kg per hectare
!
               dnileach(i) = dnileach(i) * 1.e+04 * 365.
!
! cropn is in units of kg/ha 
!
               ypnuptake = 0.0
                  do 335 j = 1, npft
                     ypnuptake = ypnuptake + cropn(i,j)
335               continue
               ypnuptake = ypnuptake + totnvegn(i) * 1.e+04

!
! calculate nitrogen balance for soil - crops - natural vegetation - inputs
! in kg ha-1
!
               snbalance(i) = taninp(i)*1.e+04 - ypnuptake - &
                              ftot(i)
            endif
100         continue
!
     return 

!	
end subroutine leaching
