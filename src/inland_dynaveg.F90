#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine dynaveg (iwest,jnorth)
! ---------------------------------------------------------------------

      use inland_parameters
      use inland_control, only: isimfire,iyear,iyrluc,nluc,isimland,iyear0
      use inland_compft
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comfire
      use inland_comcrop, only:isimagro

      implicit none

! local variables
      integer i,   &    ! gridcell counter
              j,   &
              nit, &    ! number of iterations
              niter, &  ! total number of iterations
              iwest,  & ! 1st lon index for subset
              jnorth    ! 1st lat index for subset! number of dimensions

      real*8 denswood, & ! kg/m**3
             saparea,  & ! in m**2
             sapspeed, & ! in mm/day
             sapvolume,& ! in m**3
             seedbio,  & ! correction term
             taufin,   & !
             trans,    & ! (2.5 mm/day)
             wood,     & ! total amount of woody biomass in gridcell
             rwork,    & ! working variable (1/total number of iterations)
             cdistinit   ! carbon lost to atmosphere by disturbance as estimated
                         ! in sumnow (used to compute correction term)  
      real*8 burn        ! this was formerly in fire subroutine

!      real*8 
!     >  aleaf(npoi,npft),   ! allocation fraction to leaves
!     >  aroot(npoi,npft),   ! allocation fraction to fine roots
!     >  awood(npoi,npft),   ! allocation fraction to wood
!     >  tauleaf(npft), ! turnover time of carbon in leaves (years)
!     >  tauroot(npft), ! turnover time of carbon in fine roots (years)
!     >  tauwood(npoi,npft), ! turnover time of carbon in wood (years)
!     >  tauwood0(npoi,npft) ! normal (unstressed) turnover time
      real*8 cbiolmin(lbeg:lend,npft) ! minimum leaf biomass used as seed.

! inland uses a small number of plant functional types:
!
!  1: tropical broadleaf evergreen tree
!  2: tropical broadleaf drought-deciduous trees
!  3: warm-temperate broadleaf evergreen tree
!  4: temperate conifer evergreen tree
!  5: temperate broadleaf cold-deciduous tree
!  6: boreal conifer evergreen tree
!  7: boreal broadleaf cold-deciduous tree
!  8: boreal conifer cold-deciduous tree
!  9: evergreen shrub
! 10: deciduous shrub
! 11: warm (c4) grass
! 12: cool (c3) grass
! 13: soybeans
! 14: maize 
! 15: spring and winter wheat
! 16: sugarcane
!
! ---------------------------------------------------------------------
! * * * specify biomass turnover parameters (years) * * *
! ---------------------------------------------------------------------
!       data tauleaf / 1.01,  & ! tropical broadleaf evergreen trees
!                      1.00,  & ! tropical broadleaf drought-deciduous trees
!                      1.00,  & ! warm-temperate broadleaf evergreen trees
!                      2.00,  & ! temperate conifer evergreen trees
!                      1.00,  & ! temperate broadleaf cold-deciduous trees
!                      2.50,  & ! boreal conifer evergreen trees
!                      1.00,  & ! boreal broadleaf cold-deciduous trees
!                      1.00,  & ! boreal conifer cold-deciduous trees
!                      1.50,  & ! evergreen shrubs
!                      1.00,  & ! deciduous shrubs
!                      1.25,  & ! warm (c4) grasses
!                      1.50,  & ! cool (c3) grasses
!                      999.0, & ! soybean
!                      999.0, & ! maize 
!                      999.0 /  ! wheat
!
!       data tauwood0 / 25.0, & ! tropical broadleaf evergreen trees
!                       25.0, & ! tropical broadleaf drought-deciduous trees
!                       25.0, & ! warm-temperate broadleaf evergreen trees
!                       50.0, & ! temperate conifer evergreen trees
!                       50.0, & ! temperate broadleaf cold-deciduous trees
!                      100.0, & ! boreal conifer evergreen trees
!                      100.0, & ! boreal broadleaf cold-deciduous trees
!                      100.0, & ! boreal conifer cold-deciduous trees
!                        5.0, & ! evergreen shrubs
!                        5.0, & ! deciduous shrubs
!                      999.0, & ! warm (c4) grasses
!                      999.0, & ! cool (c3) grasses
!                      999.0, & ! soybean
!                      999.0, & ! maize 
!                      999.0 /  ! wheat 

! FIXME: useless initialization (it's reset at the start of the do 100 below)
      wood = 0.001

! iteration
      niter = 10
      rwork = 1. / float(niter)

! set fire disturbance depending on isimfire value
! 0: no fire disturbance (0%/yr), 1: natural const (0.5%/yr), 2: CTEM, 3: IBIS (default 1, ignored if isimveg=0)
      if (isimfire.eq.0) then ! no fire disturbance (0%/yr)
         disturbf(:) = 0.0
      else if (isimfire.eq.1) then ! natural const (0.5%/yr)
         disturbf(:) = 0.005
      else if (isimfire.eq.2) then ! CTEM method
         disturbf(:) = pfireyr(:)
      else if (isimfire.eq.3) then ! former IBIS subroutine
         do i = lbeg, lend
            burn = firefac(i) * min (dble(1.0), totlit(i) / 0.200)
            disturbf(i) = 1.0 - exp(-0.5 * burn)
            disturbf(i) = max (dble(0.0), min (dble(1.0), disturbf(i)))
         end do
      endif

! initialize wood for gridcell
    if (isimagro .eq. 0) then

#ifndef SINGLE_POINT_MODEL
! disturbances caused by land use are set to 0, unless
! there is available data to be assimilated
      disturbl = 0.0
      if (iyear.ge.iyrluc.and.iyear.le.iyrluc+nluc.and.isimland.gt.0) then
         call rdlndtrans(iwest,jnorth)
      endif
#endif

! begin global grid
      do 100 i = lbeg, lend

! initialize wood for gridcell
         wood = 0.001
! ---------------------------------------------------------------------
! * * * initialize vegetation dynamics pools * * *
! ---------------------------------------------------------------------
!
! zero out litter fall fields
         falll(i) = 0.0
         fallr(i) = 0.0
         fallw(i) = 0.0
         caccount(i) = 0.0

! zero out carbon lost due to disturbance
         cdisturb(i) = 0.0
         cdistinit = 0.0

! set fixed disturbance regime other than fire
         disturbo(i) = 0.005

!gabriel apagar
if (i.eq.1) then
write(*,*) "falll 1   :",falll(1)
end if

#ifndef SINGLE_POINT_MODEL
! ---------------------------------------------------------------------
! * * * Balance between disturbances * * *
! ---------------------------------------------------------------------
!
      if ((disturbo(i)+disturbf(i)) .gt. 1.0) then
         disturbf(i) = 1.0 - disturbo(i)
         disturbl(i) = 0.0
      elseif ((disturbo(i)+disturbf(i)+disturbl(i)) .gt. 1.0) then
         disturbl(i) = 1.0 - disturbo(i) - disturbf(i)
      endif
#endif

! ---------------------------------------------------------------------
! * * * update npp, and pool losses  * * *
! ---------------------------------------------------------------------

! go through all the pfts
         do 110 j = 1, 12  !changed npft to 12 (number of pft) Pousa, 2013

! maintain minimum value of leaf carbon in areas where plants exist
            cbiolmin(i,j) = exist(i,j)*xminlai/specla(i,j)	! Castanho HP, 2013
!
! existence arrays applied to npp in sumnow
!        aynpp(i,j)  = exist(i,j) * aynpp(i,j)              ! Castanho HP, 2013 

! determine above-ground npp for each plant type

            ayanpp(i,j) = (aleaf(i,j) + awood(i,j)) * aynpp(i,j)    ! Castanho HP, 2013 

! determine turnover rates for woody biomass:
! if pft can exist,    then tauwood = tauwood0 (normal turnover),
! if pft cannot exist, then tauwood = taufin years (to kill off trees)
!           taufin     = 5.0
            taufin     = tauwood0(i,j) * 0.5			! Castanho HP, 2013
            tauwood(i,j) = tauwood0(i,j)-(tauwood0(i,j)-taufin) * (1.0-exist(i,j))  ! Castanho HP, 2013

! assume a constant fine root turnover time
!          tauroot(j) = 1.0

! calculate carbon lost to atmosphere by disturbance (non iterated) :
! corresponds to value calculated by sumnow and used in the instantaneous nee
! used to balance carbon

#ifndef SINGLE_POINT_MODEL
            if (isimland .eq. 1 ) then
               if (j.le.8) then
                  cdistinit = cdistinit + ((cbiol(i,j) - cbiolmin(i,j)) * &
                              (disturbf(i) + disturbl(i) + disturbo(i)) + cbiow(i,j) *  &
                              (disturbf(i) + disturbl(i) + disturbo(i)) + cbior(i,j) *  &
                              (disturbf(i) + disturbl(i) + disturbo(i)))
               else 
	                   cdistinit = cdistinit + ((cbiol(i,j) - cbiolmin(i,j)) * &
                           (disturbf(i) + disturbo(i)) + cbiow(i,j) *  &
                           (disturbf(i) + disturbo(i)) + cbior(i,j) *  &
                           (disturbf(i) + disturbo(i)))
	       endif		  
            else
               cdistinit = cdistinit + ((cbiol(i,j) - cbiolmin(i,j)) * &
                           (disturbf(i) + disturbo(i)) + cbiow(i,j) *  &
                           (disturbf(i) + disturbo(i)) + cbior(i,j) *  &
                           (disturbf(i) + disturbo(i)))
            endif
#else
            cdistinit = cdistinit + (                                     &
                           (cbiol(i,j) - cbiolmin(i,j)) * (               &
                              disturbf(i) + disturbo(i)                   &
                           ) + cbiow(i,j) * (disturbf(i) + disturbo(i)) + &
                           cbior(i,j) * (disturbf(i) + disturbo(i)))
#endif

! iteration loop
            do 10 nit = 1, niter

! determine litter fall rates
if ((nit.eq.1 .and. i.eq.1) .and. j .eq. 11) then
write(*,*) "falll antes:     ",falll(1)
write(*,*) "cbiol      :     ",cbiol(1,11)
write(*,*) "aycbiol    :     ",aycbiol(1,11)
write(*,*) "cbiolmin   :     ",cbiolmin(1,11)
write(*,*) "tauleaf    :     ",tauleaf(11)
write(*,*) "resultado  :     ",(cbiol(i,j) - cbiolmin(i,j)) / tauleaf(j)
write(*,*) "resultado2 :     ",(aycbiol(i,j) - cbiolmin(i,j)) / tauleaf(j)
write(*,*) "consumo    :     ",(numherb * 0.0001 * consherb/(1-grazef))
write(*,*) "rwork      :     ",rwork
if ((iyear - iyear0) .ge. nyrstherb) then
write(*,*) "EATING ------ EATING ------ EATING ------ EATING ------ EATING ------ EATING ------ EATING",iyear,iyear0
end if
end if


!gabriel abrahao: if isenescen is on, use yearly averages of biomass pools to calculate litter fall rates
               if (isenescen.ge.1) then
!gabriel abrahao: in case herbivory is happening, add (1-grazef) of the consumption to falll, but never add more than the current cbiol
                  if (((iyear - iyear0) .ge. nyrstherb) .and. j.eq.herbpft) then
                     falll(i) = falll(i) + min(aycbiol(i,j),((aycbiol(i,j) - cbiolmin(i,j))/tauleaf(j) + (numherb * 0.0001 * consherb/(1-grazef)))) * rwork
                  else
                     falll(i) = falll(i) + (aycbiol(i,j) - cbiolmin(i,j)) / tauleaf(j) * rwork
                  end if
                  fallr(i) = fallr(i) + aycbior(i,j) / tauroot(j) * rwork
                  fallw(i) = fallw(i) + aycbiow(i,j) / tauwood(i,j) * rwork
               else !isenesen
!gabriel abrahao: if isenescen is off, just use the current biomass pools
                  falll(i) = falll(i) + (cbiol(i,j) - cbiolmin(i,j)) / tauleaf(j) * rwork
                  fallr(i) = fallr(i) + cbior(i,j) / tauroot(j) * rwork
                  fallw(i) = fallw(i) + cbiow(i,j) / tauwood(i,j) * rwork			! Castanho HP, 2013
               end if !isenescen



if ((nit.eq.niter .and. i.eq.1) .and. j .eq. 11) then
write(*,*) "falll depois:     ",falll(1)
end if

! ---------------------------------------------------------------------
! * * * apply disturbances * * *
! ---------------------------------------------------------------------
!
! calculate biomass (vegetations) carbon lost to atmosphere   
! used to balance net ecosystem exchange

#ifndef SINGLE_POINT_MODEL
               if (isimland .eq. 1 ) then
                  if (j.le.8) then
                     cdisturb(i) = cdisturb(i) + ((cbiol(i,j) - cbiolmin(i,j)) * (disturbf(i) + &
                                   disturbl(i) + disturbo(i)) + cbiow(i,j) * &
                                   (disturbf(i) + disturbl(i) + disturbo(i)) +    &
                                   cbior(i,j) * (disturbf(i) + disturbl(i) + disturbo(i))) * rwork               

                  else
                     cdisturb(i) = cdisturb(i) + ((cbiol(i,j) - cbiolmin(i,j)) * (disturbf(i) + &
                                   disturbo(i)) + cbiow(i,j) * (disturbf(i) + disturbo(i)) +    &
                                   cbior(i,j) * (disturbf(i) + disturbo(i))) * rwork             
                  endif
               else
                  cdisturb(i) = cdisturb(i) + ((cbiol(i,j) - cbiolmin(i,j)) * (disturbf(i) + &
                                disturbo(i)) + cbiow(i,j) * (disturbf(i) + disturbo(i)) +    &
                                cbior(i,j) * (disturbf(i) + disturbo(i))) * rwork               
               endif
#else
               cdisturb(i) = cdisturb(i) + (                                  &
                                (cbiol(i,j) - cbiolmin(i,j)) * (disturbf(i) + &
                                disturbo(i)) + cbiow(i,j) * (disturbf(i) +    &
                                disturbo(i)) + cbior(i,j) * (disturbf(i) +    &
                                disturbo(i))                                  &
                             ) * rwork           
#endif

!
! ---------------------------------------------------------------------
! * * * update biomass pools  * * *
! ---------------------------------------------------------------------
!
! update carbon reservoirs using an analytical solution
! to the original carbon balance differential equation

#ifndef SINGLE_POINT_MODEL
!gabriel abrahao: isenescen only works if isimland is 0 TODO FIXME: Currently causes double accumulation if isimland.eq.1, remove accumulation from here if isenescen.ge.1
               if (isimland .eq. 1 ) then
                  if (j.le.8) then
                     cbiol(i,j) = cbiol(i,j) + (                                     &
                                  aleaf(i,j) * max (dble(0.), aynpp(i,j)) - (cbiol(i,j) - &
                                  cbiolmin(i,j)) / tauleaf(j) - (disturbf(i) +    &
                                  disturbl(i) + disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j))) * &
                                  rwork 

                     cbiow(i,j) = cbiow(i,j) + (                                    &
                                  awood(i,j) * max (dble(0.), aynpp(i,j)) - cbiow(i,j) / &
                                  tauwood(i,j) - (disturbf(i)+disturbl(i) + disturbo(i)) *       &
                                  cbiow(i,j)) * rwork

                     cbior(i,j) = cbior(i,j) + (                                    &
                                  aroot(i,j) * max (dble(0.), aynpp(i,j)) - cbior(i,j) / &
                                  tauroot(j) - (disturbf(i)+disturbl(i) + disturbo(i)) *       &
                                  cbior(i,j)) * rwork
                  else
                     cbiol(i,j) = cbiol(i,j) + (                                     &
                                  aleaf(i,j) * max (dble(0.), aynpp(i,j)) - (cbiol(i,j) - &
                                  cbiolmin(i,j)) / tauleaf(j) - (disturbf(i) +    &
                                  disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j))) * rwork 

                     cbiow(i,j) = cbiow(i,j) + (                                    &
                                  awood(i,j) * max (dble(0.), aynpp(i,j)) - cbiow(i,j) / &
                                  tauwood(i,j) - (disturbf(i)+disturbo(i)) *       &
                                  cbiow(i,j)) * rwork

                     cbior(i,j) = cbior(i,j) + (                                    &
                                  aroot(i,j) * max (dble(0.), aynpp(i,j)) - cbior(i,j) / &
                                  tauroot(j) - (disturbf(i)+disturbo(i)) *       &
                                  cbior(i,j)) * rwork 
                  endif
	       else   
!Gabriel Abrahao: Accounts for herbivore leaf consumption in ptfs .ge.9 at the end of the year if isenescen.eq.0
!gabriel apagar
if (isenescen .eq. 0) then
               if (((iyear - iyear0) .ge. nyrstherb) .and. j.eq.herbpft) then
                  cbiol(i,j) = cbiol(i,j) + (                                     &
                	       aleaf(i,j) * max (dble(0.0), aynpp(i,j)) - (cbiol(i,j) - &
                	       cbiolmin(i,j)) / tauleaf(j) - (disturbf(i) +    &
                	       disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j)) - (numherb * 0.0001 * consherb/grazef)) * rwork ! The 0.0001 converts consumption per ha to per m^2
               else
                  cbiol(i,j) = cbiol(i,j) + (                                     &
                	       aleaf(i,j) * max (dble(0.0), aynpp(i,j)) - (cbiol(i,j) - &
                	       cbiolmin(i,j)) / tauleaf(j) - (disturbf(i) +    &
                	       disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j))) * rwork ! The 0.0001 converts consumption per ha to per m^2
               end if

                  cbiow(i,j) = cbiow(i,j) + (				      &
                       	       awood(i,j) * max (dble(0.0), aynpp(i,j)) - cbiow(i,j) / &
                	       tauwood(i,j) - (disturbf(i)+disturbo(i)) * cbiow(i,j)) * rwork

                  cbior(i,j) = cbior(i,j) + (				      &
                    	       aroot(i,j) * max (dble(0.0), aynpp(i,j)) - cbior(i,j) / &
                	       tauroot(j) - (disturbf(i)+disturbo(i)) * cbior(i,j)) * rwork
!gabriel abrahao: If isenescen is 1, the C balance equation is solved elsewhere (currently in main_offline) and only the disturbances should be done here
else
                  cbiol(i,j) = cbiol(i,j) + (  0.0                                   &
                	        - (disturbf(i) +    &
                	       disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j))) * rwork 

                  cbiow(i,j) = cbiow(i,j) + ( 0.0				      &
                	        - (disturbf(i)+disturbo(i)) * cbiow(i,j)) * rwork

                  cbior(i,j) = cbior(i,j) + ( 0.0				      &
                	        - (disturbf(i)+disturbo(i)) * cbior(i,j)) * rwork
end if
!gabriel apagar
if (i.eq.4 .and. j.eq.11 .and. nit.eq.10) then
	write(*,*) "WARNING: Changes applied to carbon balance in dynaveg"
	write(*,*) "cbiol:           ",cbiol(i,j)
	write(*,*) "cbiow:           ",cbiow(i,j)
	write(*,*) "cbior:           ",cbior(i,j)
	!write(*,*) "deltacbiol: ",aleaf(i,j) * max (dble(0.), aynpp(i,j)) - (cbiol(i,j) -cbiolmin(i,j)) / tauleaf(j) - (disturbf(i) +disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j))  
	write(*,*) "aynpp:           ",aynpp(i,j)
	write(*,*) "cbiol/tauleaf:   ",(cbiol(i,j) - cbiolmin(i,j))/tauleaf(j)
	write(*,*) "disturbf:        ",disturbf(i)
	write(*,*) "disturbo:        ",disturbo(i)
	!write(*,*) "cbiolmin:   ",cbiolmin(i,j)

!open(1991,file="/home/gabriel/doutorado/inlands/testes/pastagem_12AGO_mother_rd3d_teste/yearly_cbiol.csv",access="APPEND")
!write(1991,*) iyear,cbiol(4,11)
end if

!!end gabriel apagar
               endif
#else
               cbiol(i,j) = cbiol(i,j) + (                                     &
                               aleaf(i,j) * max (dble(0.), aynpp(i,j)) - (cbiol(i,j) - &	! Castanho HP, 2013
                               cbiolmin(i,j)) / tauleaf(j) - (disturbf(i) +    &
                               disturbo(i)) * (cbiol(i,j) - cbiolmin(i,j))     &
                            ) * rwork 

               cbiow(i,j) = cbiow(i,j) + (                                    &
                               awood(i,j) * max (dble(0.), aynpp(i,j)) - cbiow(i,j) / &		! Castanho HP, 2013
                               tauwood(i,j) - (disturbf(i)+disturbo(i)) *       &			! Castanho HP, 2013
                               cbiow(i,j) &
                            ) * rwork

               cbior(i,j) = cbior(i,j) + (                                    &
                               aroot(i,j) * max (dble(0.), aynpp(i,j)) - cbior(i,j) / &		! Castanho HP, 2013
                               tauroot(j) - (disturbf(i)+disturbo(i)) *       &
                               cbior(i,j) &
                       ) * rwork
#endif

10			continue ! end of iteration loop

!gabriel abrahao:FIXME: this may be a dangerous approach. If isenescen is 1, biomass pools values are not as meaningful in an yearly context because in this case we are looking at the december value instead of an actual yearly value. So we save the pools to dummy variables here and use the average yearly variables for the rest of dynaveg, recovering the original monthly values in the end
            if (isenescen.ge.1) then
               dumcbiol(i,j) = cbiol(i,j)
               dumcbior(i,j) = cbior(i,j)
               dumcbiow(i,j) = cbiow(i,j)
               dumplai(i,j) = plai(i,j)

               cbiol(i,j) = aycbiol(i,j)
               cbior(i,j) = aycbior(i,j)
               cbiow(i,j) = aycbiow(i,j)
               plai(i,j)  = ayplai(i,j)
            end if

            if (j.le.8) wood = wood + max (dble(0.0), cbiow(i,j))
               seedbio = max(dble(0.),(cbiolmin(i,j) - cbiol(i,j)))

! account for negative biomass in nee (caccount > 0: carbon that has been accounted for 
! as absorbed in the fluxes but that is not accounted for in the calculation of biomass 
! pools ==> has to be released to atmosphere)
            caccount(i) = caccount(i)+seedbio-min(dble(0.0), cbiow(i,j))-min(dble(0.0), cbior(i,j))

! constrain biomass fields to be positive
            cbiol(i,j) = max (cbiolmin(i,j), cbiol(i,j))
            cbiow(i,j) = max (dble(0.0), cbiow(i,j))
            cbior(i,j) = max (dble(0.0), cbior(i,j))

! update vegetation's physical characteristics
            biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
            if (isenescen .eq. 0) then
               plai(i,j)    = cbiol(i,j) * specla(i,j)			! Castanho HP, 2013
            else
               !gabriel.abrahao@ufv.br: If isenescen is not zero, this should be done somewhere else
            end if
110      continue

! ---------------------------------------------------------------------
! * * * update annual npp, lai, and biomass * * *
! ---------------------------------------------------------------------

! Disturbance can't result in negative biomass. caccount account for the 
! carbon not to be removed by the disturbance.
!       cdisturb(i) = cdisturb(i) - caccount(i)
!
! determine total ecosystem positive npp (changed by exist at begin
! of subroutine). Different from sum of monthly and daily npp)
         aynpptot(i) = max(dble(0.0),aynpp(i,1))  + max(dble(0.0),aynpp(i,2))  + &
                       max(dble(0.0),aynpp(i,3))  + max(dble(0.0),aynpp(i,4))  + &
                       max(dble(0.0),aynpp(i,5))  + max(dble(0.0),aynpp(i,6))  + &
                       max(dble(0.0),aynpp(i,7))  + max(dble(0.0),aynpp(i,8))  + &
                       max(dble(0.0),aynpp(i,9))  + max(dble(0.0),aynpp(i,10)) + &
                       max(dble(0.0),aynpp(i,11)) + max(dble(0.0),aynpp(i,12))

! adjust annual net ecosystem exchange (calculated in stats.f) 
! by new value of npp (depending on exist), andloss of carbon to
! atmosphere due to biomass burning (fire)
!       ayneetot(i) = aynpptot(i) - ayco2mic(i) - cdisturb(i) + caccount(i)

! the fact that only >0 npp was used has also to be accounted for (only 
! possible 1ce a year) difference between disturbance calculated in dynaveg and
! sumnow has also to be accounted for
         caccount(i) = caccount(i) -                                 &
                       min(dble(0.0),aynpp(i,1))  - min(dble(0.0),aynpp(i,2))  - &
                       min(dble(0.0),aynpp(i,3))  - min(dble(0.0),aynpp(i,4))  - &
                       min(dble(0.0),aynpp(i,5))  - min(dble(0.0),aynpp(i,6))  - &
                       min(dble(0.0),aynpp(i,7))  - min(dble(0.0),aynpp(i,8))  - &
                       min(dble(0.0),aynpp(i,9))  - min(dble(0.0),aynpp(i,10)) - &
                       min(dble(0.0),aynpp(i,11)) - min(dble(0.0),aynpp(i,12)) - &
                       (cdisturb(i) - cdistinit)

! adjust annual net ecosystem exchange (calculated in sumyear)
! by corrections computed in dynaveg
         ayneetot(i) = ayneetot(i) + caccount(i)

! determine total ecosystem above-ground npp
         ayanpptot(i) = ayanpp(i,1)  + ayanpp(i,2)  + ayanpp(i,3) + &
                        ayanpp(i,4)  + ayanpp(i,5)  + ayanpp(i,6) + &
                        ayanpp(i,7)  + ayanpp(i,8)  + ayanpp(i,9) + &
                        ayanpp(i,10) + ayanpp(i,11) + ayanpp(i,12)

! update total canopy leaf area
         totlaiu(i) = plai(i,1) + plai(i,2) + plai(i,3) + plai(i,4) + &
                      plai(i,5) + plai(i,6) + plai(i,7) + plai(i,8)
         totlail(i) = plai(i,9) + plai(i,10) + plai(i,11) + plai(i,12)

! update total biomass
         totbiou(i) = biomass(i,1) + biomass(i,2) + biomass(i,3) + &
                      biomass(i,4) + biomass(i,5) + biomass(i,6) + &
                      biomass(i,7) + biomass(i,8)
         totbiol(i) = biomass(i,9)  + biomass(i,10) + biomass(i,11) + &
                      biomass(i,12)

! ---------------------------------------------------------------------
! * * * update fractional cover and vegetation height parameters * * *
! ---------------------------------------------------------------------

! update fractional cover of forest and herbaceous canopies:
         fu(i) = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
         fl(i) = totlail(i) / 1.0

! apply disturbances to fractional cover

#ifndef SINGLE_POINT_MODEL
         if (isimland .eq. 1 ) then
            fu(i) = fu(i) * (1. - disturbf(i) - disturbl(i) + disturbo(i))
            fl(i) = fl(i) * (1. - disturbf(i) - disturbl(i) + disturbo(i))
         else
            fu(i) = fu(i) * (1. - disturbf(i) - disturbo(i))
            fl(i) = fl(i) * (1. - disturbf(i) - disturbo(i))	 
	 endif
#else
         fu(i) = fu(i) * (1. - disturbf(i) - disturbo(i))
         fl(i) = fl(i) * (1. - disturbf(i) - disturbo(i))
#endif

! constrain the fractional cover
         fu(i) = max(dble(0.25), min(dble(0.975), fu(i)))
         fl(i) = max(dble(0.25), min(dble(0.975), fl(i)))

! annual update upper canopy height parameters
! should be calculated based on vegetative fraction and not the
! average over the entire grid cell
         zbot(i,2) = 3.0
         ztop(i,2) = max(zbot(i,2) + 1.00, 2.50 * totbiou(i) / fu(i) * 0.75)

! ---------------------------------------------------------------------
! * * * update stem area index and sapwood fraction * * *
! ---------------------------------------------------------------------

! estimate stem area index (sai) as a fraction of the lai
         sai(i,1) = 0.050 * totlail(i)
         sai(i,2) = 0.250 * totlaiu(i)

! estimate sapwood fraction of woody biomass
         sapspeed  = 25.0                        ! (m/day)
         trans     = 0.0025                      ! (2.5 mm/day) 
         saparea   = (trans / sapspeed)          ! m**2
         sapvolume = saparea * ztop(i,2) * 0.75  ! m**3
         denswood  = 400.0                       ! kg/m**3
         sapfrac(i) = min (dble(0.50), max (dble(0.05), sapvolume * denswood / wood))
!gabriel.abrahao@ufv.br apagar
!totbiol(i) = cbiol(i,11)
!write(*,*) "WARNING: This version sets biomass to be only cbiol, see gabriel.abrahao@ufv.br in inland_dynaveg"
if (i.eq.4) then
! write(*,*) 
! write(*,*) "plai:		",plai(i,11)
! write(*,*) "cbiol:		",cbiol(i,11)
! write(*,*) "biomass:		",biomass(i,11)
! write(*,*) "totbiol:		",totbiol(i)
! write(*,*) "ciobl/biomass:	",cbiol(i,11) / biomass(i,11)
! write(*,*) "plai/cbiol:	",plai(i,11) / cbiol(i,11)
!gabriel apagar
	write(*,*) "YEARLY -- YEARLY --YEARLY --YEARLY --YEARLY --YEARLY --YEARLY --YEARLY --YEARLY --YEARLY --YEARLY --YEARLY"
	write(*,*) "ayco2mic:  ", ayco2mic(4)
	write(*,*) "cdisturb:  ", cdisturb(4)
	write(*,*) "caccount:  ", caccount(4)
end if

100   continue
       endif  ! check for crop existence 
! ---------------------------------------------------------------------
! * * * map out vegetation classes for this year * * *
! ---------------------------------------------------------------------
      call vegmap
!gabriel abrahao:FIXME: this may be a dangerous approach. If isenescen is 1, biomass pools values are not as meaningful in an yearly context because in this case we are looking at the december value instead of an actual yearly value. So the pools were saved after applying the disturbances to dummy variables and we used the average yearly variables for the rest of dynaveg, recovering the original monthly values here
            if (isenescen.ge.1) then
               cbiol(:,:) = dumcbiol(:,:)
               cbior(:,:) = dumcbior(:,:)
               cbiow(:,:) = dumcbiow(:,:)
               plai(:,:)  = dumplai(:,:) 
            end if


      return
end  subroutine dynaveg
