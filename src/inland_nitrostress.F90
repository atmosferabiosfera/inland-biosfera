#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine nitrostress(istep,iday) 
! ---------------------------------------------------------------------
!
! this subroutine calculates the effects of the amount of available 
! inorganic nitrogen on carbon assimilation in crops 
!
! strictly speaking, stressn* is multiplied to the vmax parameters 
! used in the photosynthesis calculations
!
! common blocks
!
!
      use inland_parameters
      use inland_control,only:imonth
      use inland_com1d
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
!
! local variables
!
      real*8 tnsupply, &    ! total nitrogen supply from soil to plant 
             gfn,      &    ! function controlling nitrogen uptake
             wsupply,  &    ! supply of water to plant from transpiration stream (kg/h20/m2/s)
!    	     alphacc,  &
             awc,      &
             sumnval,  &
             fnmin,    &
             wc,       &
             dmplant,  &
             fdmleaf,  &
             fdmstem,  &
             fdmroot,  &
             fdmgrain, &
             dmres,    &
             flres,    &
             fsres,    &
             tngrain,  &
             tnleaf,   &
             tnstem,   &
             tnroot,   &
             fnsmax,   &
             fnrmax,   &
             fnpmax,   &
             f1,f2,f3
!
! biological fixation
!
      real*8 gs,     &
             fxg,    &
             rdepth, &
             fc,     &
             wp,     &
             sm,     &
             fxw,    &
             rd,     &
             fxn,    &
             fxr,    &
             fnresidue
!
      real*8 plantn(lbeg:lend)    ! plant available nitrogen total
!
!     real*8 fnlfmx(npft), &
!    	     fngrmx(npft), &
!    	     sratio(npft), &
!    	     rratio(npft), &
!    	     fnopt(npft)

      integer iday,   &
              istep,  &
              iter,   &
              iter2,  &
              i,j,k
!
!    cfrac(13)    = 0.50  ! fraction of dry matter that is carbon   
!    cfrac(14)    = 0.50    
!    cfrac(15)    = 0.45    
!
!    fnlfmx(13)   = 0.025 ! average leaf nitrogen content  
!    fnlfmx(14)   = 0.013 ! average leaf nitrogen content
!
!    if (iwheattype .eq. 1) fnlfmx(15)   = 0.008 ! average leaf nitrogen content  
!    if (iwheattype .eq. 2) fnlfmx(15)   = 0.010 ! average leaf nitrogen content 
!
     if ( iwheattype .ne. 0 ) then
          fnlfmx(15) = fnlfmxw(iwheattype)
!
!    fnlfmx(15)   = 0.060 ! average leaf nitrogen content  
!
!    fngrmx(13)   = 0.040 ! maximum allowed grain nitrogen fraction
!    fngrmx(14)   = 0.017
!
!    if (iwheattype .eq. 1) fngrmx(15)   = 0.0300 ! maximum allowed grain nitrogen fraction
!    if (iwheattype .eq. 2) fngrmx(15)   = 0.0225 ! maximum allowed grain nitrogen fraction
!
          fngrmx(15) = fngrmxw(iwheattype)   
!
!    sratio(13)   = 0.40  ! typical ratio of stem to leaf nitrogen content
!    sratio(14)   = 0.05  
!    sratio(15)   = 0.40  ! typical ratio of stem to leaf nitrogen content
!
!    rratio(13)   = 0.75  ! typical ratio of root to leaf nitrogen content
!    rratio(14)   = 0.75  ! typical ratio of root to leaf nitrogen content
!    rratio(15)   = 1.00  ! typical ratio of root to leaf nitrogen content
!
!    fnopt(13) = 0.0075   ! minimum leaf nitrogen concentration - stress onset - soybean
!    fnopt(14) = 0.03125  ! leaf nitrogen concentration - stress onset - maize 
!    fnopt(14) = 0.02850  ! leaf nitrogen concentration - stress onset - maize 
!
!    if (iwheattype .eq. 1) fnopt(15) = 0.0175   ! stress onset for leaf nitrogen - spring wheat 
!    if (iwheattype .eq. 2) fnopt(15) = 0.0100   ! stress onset for leaf nitrogen - winter wheat 
!
          fnopt(15) = fnoptw(iwheattype)
     else
          fnlfmx(15) = 0
          fngrmx(15) = 0     
          fnopt(15) =  0
     endif 
! 
     do 100 i = lbeg, lend
!
!      alphac = 0.0005   ! minimum uptake rate applied to nitrogen (mm/timestep) 
!      gnmin  = 0.001	 ! minimum nitrogen fraction allowed in grain  
!      smax   = 1.05	 ! maximum nitrogen stress factor value     
!      availn = 1.0	 ! scaling variable to adjust for plant capability in capturing
			  ! nitrogen from pools - above 1.0 means plant can take
			  ! some up in excess of what is in the transpiration stream 
	fngrain(i,13) = 0.035  
	fngrain(i,14) = 0.013
	fngrain(i,15) = 0.020  
	fngrain(i,16) = 0.013
!
	tnsupply  = 0.0
	awc	  = 0.0
	sumnval   = 0.0
!
! initialize layer nitrogen uptake
!
	do 180 k = 1, nsoilay
!
	   tnuptake(i,k) = 0.0
	   anuptake(i,k) = 0.0
!
180	continue
!
        do 200 j = scpft, ecpft 
!
           if (exist(i,j) .eq. 1.0) then
	      if (croplive(i,j) .eq. 1.0) then 
!                 cnmax       = 95 
                 fnmin       = cfrac(j) / cnmax
                 stressn(i,j)= 1.0
                 gfn         = 1.0
!
! calculate the total nitrogen supply rate (kg/m2/day) for each soil layer based on
! 1) the total daily water uptake for crops (lower canopy) (mm/day - stats.f)
! 2) total available nitrogen pool to roots (soil solution) (kg-no3/m2)
! 3) and available water content (mm)
!
! NOTE:  at this time, logic in IBIS cannot be used to determine
! the uptake of nitrogen for each specific pft (mixed in each grid
! cell because upsoil is for the entire lower canopy...will have
! to weight it for now on the lai of that pft [frac(i,j)] 
!
! since it is being called each timestep, use instantaneous values
! of soil ice and moisture
!  
	   do 210 k = 1, nsoilay 
!
! calculate water content in each layer - based on EPIC parameterizations
! that look at actual water content in mm and not available water 
!
	      wc = max(0.0, (wisoi(i,k) +  &
		   (1.0 - wisoi(i,k)) *    &
		   wsoi(i,k))) * hsoi(k) * & 
		   poros(i,k) * 1000
!
!
! alphac is minimum uptake rate to account for nitrogen usage
! even when transpiration is small (plant still able to take up
! nitrogen)
!
! allow plant to take up nitrogen in excess of the
! transpiration stream at low rates early in the
! season 
!
! supply of nitrogen to crops is from roots corresponding to [l]ower
! canopy in model - since this routine is being called each timestep
! use the value of upsoil(i,k) from canopy.f rather than the value from
! stats.f which is the daily average (adupsoil)   
! upsoil is in units of mm/m2/s of transpiration
!
	      wsupply = upsoil(i,k) * dtime
!
! make sure that water content is not zero - 
! set to small limit
!
	      wc = max(1.0, wc)
!
! the total nitrogen uptake from the layer comes from the total n
! in the layer both in solution and in soil - leachable n is only
! that portion that is in the solution
!
! value of tnuptake for this layer is used in subroutine in solute
! leaching algorithm as a net sink of nitrogen to the layer 
!
! only allow uptake in layers that have roots
! make sure nitrogen uptake only occurs while crop is active
! the minimum rate will only be applied when the plant is not
! experiencing moisture stress
!
	      if (froot(k,1) .gt. 0.005 .and. tnpptot(i) .gt. 0.0) then
		 tnuptake(i,k) = max(alphac * stressl(i,k), wsupply) * availn * &
				(smsoil(i,k) + smsoln(i,k)) 

	      else
		 tnuptake(i,k) = 0.0
	      endif
!
210	   continue
!
!	   endif
	   
	   if (aybprod(i,j) .gt. 0.0  .and. &
	       aylprod(i,j) .gt. 0.0) then
!
! for the purpose of dealing with total nitrogen uptake,
! we have to use year to date total carbon production
! in these equations because some root and leaf biomass
! has been adjusted due to phenology in the model
!
	      dmplant	=  aybprod(i,j)  / cfrac(j)
	      fdmleaf	=  (aylprod(i,j) / cfrac(j)) / dmplant
	      fdmstem	=  (cbios(i,j)   / cfrac(j)) / dmplant
	      fdmroot	=  (ayrprod(i,j) / cfrac(j)) / dmplant
	      fdmgrain  =  (cbiog(i,j)   / cfrac(j)) / dmplant
!
	      dmres	=  (aylprod(i,j) + cbios(i,j)) / cfrac(j) 
	      flres	=  (aylprod(i,j) / cfrac(j)) / dmres
	      fsres	=  (cbios(i,j)   / cfrac(j)) / dmres
!
	      fnplant(i,j)   =  max(0.0, totnuptake(i,j) / dmplant)

!
! maintain minimum nitrogen concentration in leaf and stem (potential residue)
!
	      iter  = 0
	      iter2 = 0
!
400           fnleaf(i,j) = (fnplant(i,j) - fngrain(i,j) * fdmgrain) / &
                            (fdmleaf + sratio(j) * fdmstem + &
                            rratio(j) * fdmroot)
!
            fnleaf(i,j) = max(0.0, fnleaf(i,j))
            fnstem(i,j) = fnleaf(i,j) * sratio(j)
            fnroot(i,j) = fnleaf(i,j) * rratio(j)
!
	    fnresidue	= fnleaf(i,j) * flres + fnstem(i,j) * fsres
	    if (fnresidue .gt. fnmin) iter2 = 1
	    if (fnresidue .le. fnmin) iter  = 1 

!	      if (iter .eq. 1 .and. fngrain(i,j) .gt. gnmin &
!			     .and. fdmgrain	.gt. 0.0 &
!			     .and. iter2	.eq. 0) then 
!		     fngrain(i,j) = max(gnmin, fngrain(i,j) * 0.99)
!		     write(*,*) 'taking from grain', fngrain(i,j)
!		     goto 400
!	  
	      if (iter2 .eq. 1 .and. fngrain(i,j) .lt. fngrmx(j) &
			       .and. fdmgrain .gt. 0.0 &
			       .and. iter .eq. 0) then 
		    fngrain(i,j) = min(fngrmx(j), fngrain(i,j) * 1.01)
		    goto 400
	   endif
!
! ------------------------------------------------------------------------------
! calculate nitrogen content in various pools
!
	   tngrain = fngrain(i,j) * fdmgrain * dmplant
	   tnleaf  = fnleaf(i,j)  * fdmleaf  * dmplant
	   tnstem  = fnstem(i,j)  * fdmstem  * dmplant
	   tnroot  = fnroot(i,j)  * fdmroot  * dmplant
!
	   tnplant(i,j) = tngrain + tnleaf + tnstem + tnroot
!
	   fnsmax   = sratio(j)  * fnlfmx(j)
	   fnrmax   = rratio(j)  * fnlfmx(j) 
!
	   fnpmax   = fnlfmx(j) * fdmleaf + fnsmax * fdmstem + &
		      fnrmax * fdmroot + fngrmx(j) * fdmgrain	    
!
!
! calculate function controlling rate of nitrogen uptake
! based on plant nitrogen concentration and maximum value
! there is a chance early in growing season that fnplant
! could be higher than fnpmax - thus maximize the below
! fraction to be .le. 1.0 
!
	   gfn = 1. - min(1.0, fnplant(i,j) / fnpmax) ** 1.0
!
! calculate the annual running total of the actual nitrogen
! uptake for all pfts in lower canopy (crops)
!
! adjust nitrogen uptake by plant by gfn factor - equal in each layer
!
	   do 310 k = 1, nsoilay 
	      anuptake(i,k) = tnuptake(i,k) * gfn
	      tnsupply      = tnsupply + anuptake(i,k)
310	   continue 
!
	   totnuptake(i,j) = totnuptake(i,j) + tnsupply 
!
	   totnuptake(i,j) = max(0.0, min((1.0 - cfrac(j)) * dmplant, totnuptake(i,j)))   
!
! calculate stress parameter using rectangular hyperbola which
! relates leaf nitrogen concentration to vmax	     
!
! rectangular hyperbola 
! f1 and f2  control the shape of the stress response function 
! which spans from 0 to 1.0 for leaf n concentrations from 0 to 4 percent 
!
! s-shaped curve nitrogen limitation effect from epic model
!
	   f1 = 8.5 
	   f2 = 11.0 
!
! ratio of leaf nitrogen concentration to the optimal maximum
! for corn/maize
	   f3 = 2 * (fnleaf(i,j) / fnopt(j))
!
	   stressn(i,j) = min (smax, (f3 / (f3 + exp(f1 - f2 * f3))))
!
	   stressn(i,j) = max (0.10, stressn(i,j))
!
!----------------------------------------------------------------------
! biological fixation of nitrogen through symbiosis in soybeans
!----------------------------------------------------------------------
!
! Key reference:
! M. Cabelguenne et al., Agricultural systems 60: 175-196, 1999.
! this module is taken from the new epicphase model
!
! the amount of daily n-fixation by the plant is based on a fraction
! of the total daily n-uptake.  It is controlled by three main factors:
!
! * growth stage of the crop (0-1) *
! * soil moisture	     (0-1) * 
! * nitrogen in rooting zone (0-1) * 
!
! the growth stage factor (fxp) inhibits fixation in young plants
! and old plants, and peaks between 30-55% of the crop cycle
!
! the soil water content factor (fxw) reduces n-fixation when the 
! water content in the top 0.3 m is less than 85% of field capacity
!
! the soil nitrogen (plant available) factor (fxn) reduces n-fixation
! when the nitrogen amount in the root zone is greater than 100 kg/ha  
!
!  
	   if (j .eq. 13) then
!
! calculate growth stage and factor (fraction of total average gdd)
!
	      gs = hui(i,j) / gddmaturity(i,j)
!
	      if     (gs .le. 0.15 .or.  gs .ge. 0.75) then
		 fxg = 0.0
	      elseif (gs .gt. 0.15 .and. gs .le. 0.30) then
		 fxg = 6.67 * gs - 1.0
	      elseif (gs .gt. 0.30 .and. gs .le. 0.55) then
		 fxg = 1.0
	      else
		 fxg = 3.75 - 5.0 * gs   
	      endif
!
! calculate effect of soil moisture in top 25-30 cm
!	     rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3))  
	      rdepth = 1. / (hsoi(1) + hsoi(2))  
!
	      fc = 0.0
	      wp = 0.0
	      sm = 0.0
	      plantn(i) = 0.0
!
	      do 220 k = 1, 2
		 fc = fc + sfield(i,k) * hsoi(k)
		 wp = wp + swilt(i,k)  * hsoi(k)
		 sm = sm + wsoi(i,k)   * hsoi(k)
!
! calculate available plant nitrogen total in these layers
!
		 plantn(i) = plantn(i) + smsoil(i,k) + smsoln(i,k)
!
220	      continue 
!
	      fc = fc * rdepth
	      wp = wp * rdepth
	      sm = sm * rdepth
	      sm = min(sm, 0.85 * (fc - wp) + wp)
!
	      fxw = (sm - wp) / (0.85 * (fc - wp)) 
!
! calculate effect of plant available nitrogen pool
! 
	      rd = 1.0   ! rooting depth in meters

	      if     (plantn(i) .gt. 0.0300) then
		 fxn = 0.0
	      elseif (plantn(i) .le. 0.0100) then
		 fxn = 1.0
	      else
		 fxn = 1.5 - 0.005 * (plantn(i) * 10000) / rd
	      endif
!
! equation for fxn has to be in kg/ha nitrogen and meters for rd 
!
! CJK replaced 2/1/2006  
! the 1.5 at the end of the following equation was to increase
! the annual n-fixation rates by 50% because of simulated errors
! and a low bias compared to observations
!
	      fxr  = min(1.0, fxw, fxn) * fxg 
!
	      fixn(i,j) = fxr * tnsupply 
!	      fixn(i,j) = 0.0 
! 
! fixn is thus calculated each timestep
!
! update plant available nitrogen pool for fixation from soybean
!
	      else  ! non-soybean crop - no fixation
!
		 fixn(i,j) = 0.0
!
	      endif    ! soybeans only fix nitrogen
!
	      totnfix(i,j) = totnfix(i,j) + fixn(i,j)
!
	      do 240 k = 1, nsoilay 
!
! critical: what quantity do we add the nitrogen fixation
! to?  previous? depends on what order subroutines are called
! 
		 smsoil(i,k) = smsoil(i,k) + fixn(i,j) * froot(k,1)
!
240	      continue
!
	   endif   ! production gt 0
	endif	 ! crop plant  
     endif     ! crop existence 
!
200  continue	 ! crop pft 
!
100  continue	 ! grid cell 
!
! return to program 
!
      return
!      
end subroutine nitrostress

