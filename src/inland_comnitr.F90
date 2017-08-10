#include "inland_config.h"
module inland_comnitr
!
      implicit none
      public
      save
!
! -------
! comnitr
! -------
! nitrogen related variables for crops, leaching and drainage
!
      real*8, dimension(:), allocatable :: aplantn,                               &
                                           totimm, totmin, totnrel, ctot, ctoti,  &
                                           ftot, tsinp, tslay, taninp, tsnimm,    &
                                           tsnmob, dtnleach, dnileach, tpnuptake, &
                                           totnvegn, drntot, ddrn, concn, assimn, & 
                                           ydeposn, yfixsoin,yno3leach,snbalance, &
                                           fixsoin
! 
!  aplantn(npoi):		available inorganic nitrogen in soil profile for plant growth (kg_n m-2 y-1) 
!  totimm(npoi):		total immobilized nitrogen in timestep (kg_n m-2 timestep-1) 
!  totmin(npoi):		total mineralized nitrogen in timestep (kg_n m-2 timestep-1) 
!  totnrel(npoi):		total mineralized/immobilized nitrogen in timestep (non-microbial) (kg_n m-2 timestep-1) 
!  ctot(npoi):		total inorganic nitrogen in soil profile (kg n  m-2)
!  ctoti(npoi):		initial total inorganic nitrogen in profile at beginning of day (kg n m-2)
!  ftot(npoi):		annual total inorganic nitrogen leached from soil profile (kg n  m-2 y-1)
!  tsinp(npoi):		daily total inorganic nitrogen inputs (kg n m-2 d-1)
!  tslay(npoi):		daily total inorganic nitrogen inputs to a soil depth determined in solute.f
!  taninp(npoi):		total annual inorganic nitrogen inputs to soil (kg n m-2 y-1)
!  tsnimm(npoi):		total soil inorganic nitrogen in immobile pool (kg n m-2)  
!  tsnmob(npoi):		total soil inorganic nitrogen in mobile pool (kg n m-2)  
!  dtnleach(npoi):		daily total inorganic nitrogen leached from entire profile (kg (nh4 + no3) m-2 d-1)
!  dnileach(npoi):		daily rate of nitrate-nitrogen leached from profile (kg no3 m-2 y-1)
!  tpnuptake(npoi):		daily total plant nitrogen uptake (kg n m-2 d-1) 
!  totnvegn(npoi):		annual total inorganic nitrogen uptake by natural vegetation (kg n m-2 y-1) 
!  drntot(npoi):		annual total drainage through the soil profile (mm y-1)
!  ddrn(npoi):		daily total drainage through specified layer in soil profile (mm/day)
!  concn(npoi):		nitrate concentration at specified depth in solute.f (mg/liter or ppm)
!  assimn(npoi):		annual total nitrogen assimilated by natural vegetation (kg n m-2 y-1)
!  ydeposn(npoi):		annual total nitrogen deposition (wet and dry) (kg n m-2 y-1)
!  yfixsoin(npoi):		annual total nitrogen fixation by natural vegetation (kg n m-2 y-1) 
!  yno3leach(npoi):		annual total nitrate-n leached from profile (kg n m-2 y-1)
!  snbalance(npoi):		annual soil nitrogen balance calculation (kg n ha-1)
!  fixsoin(npoi):		general nitrogen fixation value from vegetation (kg n m-2 day-1)
!

      real*8, dimension(:,:), allocatable :: totnuptake, stressn, totnfix,     &
                                             fixn, fertnitro

!  totnuptake(npoi,npft): annual total nitrogen uptake (kg_n m-2 y-1) 
!  stressn(npoi,npft):    stress factor applied to vmax based on leaf nitrogen content in crops (dimensionless)
!  totnfix(npoi,npft):    annual total nitrogen fixation through symbiosis (kg_n m-2 y-1) 
!  fixn(npoi,npft): 	  timstep total nitrogen fixation 
!  fertnitro(npoi,npft):  nitrogen added as fertilizer to each crop

      real*8, dimension(:,:), allocatable :: daydrn, daynconc

!  daydrn(npoi,366): 	daily drainage at specified depth in solute.f (mm/day)
!  daynconc(npoi,366):  daily average nitrate concentration at specified depth in solute.f (mg/liter)  
!
end module inland_comnitr
