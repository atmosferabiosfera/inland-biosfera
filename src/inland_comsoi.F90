#include "inland_config.h"
module inland_comsoi

      implicit none
      public
      save
!
! ------
! comsoi
! ------
!
      real*8 wpudmax, zwpmax, bperm, h20
!
!     wpudmax                  ! normalization constant for puddles (kg m-2)
!     zwpmax                   ! assumed maximum fraction of soil surface 
!                              ! covered by puddles (dimensionless)
!     bperm                    ! lower boundary condition for the soil (permeability)
!                              ! bperm = 0.00  ! e.g. fully impermeable base
!                              ! bperm = 1.00  ! e.g. fully permeable base
!     h20                      ! C leaching fudge factor
!
      real*8, dimension(:), allocatable :: wpud, wipud, z0soi, albsav, albsan, stresstl, stresstu, &
                                           heati, heatg, hvasug, hvasui, soihfl, tg, ti
					   
!
!     wpud(npoi)               ! liquid content of puddles per soil area (kg m-2)
!     wipud(npoi)              ! ice content of puddles per soil area (kg m-2)
!     z0soi(npoi)              ! roughness length of soil surface (m)
!     albsav(npoi)             ! saturated soil surface albedo (visible waveband)
!     albsan(npoi)             ! saturated soil surface albedo (near-ir waveband)
!     stresstl(npoi)           ! sum of stressl over all 6 soil layers (dimensionless)
!     stresstu(npoi)           ! sum of stressu over all 6 soil layers (dimensionless)
!     heati(npoi)              ! net heat flux into snow surface (W m-2)
!     heatg(npoi)              ! net heat flux into soil surface (W m-2)
!     hvasug(npoi)             ! latent heat of vap/subl, for soil surface (J kg-1)
!     hvasui(npoi)             ! latent heat of vap/subl, for snow surface (J kg-1)
!     soihfl(npoi)             ! soil heat flux (W m-2)
!     tg(npoi)                 ! soil skin temperature (K)
!     ti(npoi)                 ! snow skin temperature (K)
!
      real*8, dimension(:), allocatable :: hsoi
!
!     hsoi(nsoilay+1)          ! soil layer thickness (m)
!
      real*8, dimension(:,:), allocatable :: tsoi, wsoi, wisoi, consoi, csoi, hydraul, suction, &
                                             bex, sfield, swilt, rhosoi, poros, porosflo, sand, &
                                             clay, stressl, stressu, upsoiu, upsoil, tnuptake,  &
                                             anuptake, smsoil, smsoln, fout, nout, csoln,       &
                                             domtext, fracsand, fracsilt,     &
                                             fracclay, swater, sice
!
!     tsoi(npoi,nsoilay)       ! soil temperature for each layer (K)
!     wsoi(npoi,nsoilay)       ! fraction of soil pore space containing liquid water
!     wisoi(npoi,nsoilay)      ! fraction of soil pore space containing ice
!     consoi(npoi,nsoilay)     ! thermal conductivity of each soil layer (W m-1 K-1)
!     csoi(npoi,nsoilay)       ! specific heat of soil, no pore spaces (J kg-1 deg-1)
!     hydraul(npoi,nsoilay)    ! saturated hydraulic conductivity (m/s)
!     suction(npoi,nsoilay)    ! saturated matric potential (m-h2o)
!     bex(npoi,nsoilay)        ! exponent "b" in soil water potential
!     sfield(npoi,nsoilay)     ! field capacity soil moisture value (fraction of pore space)
!     swilt(npoi,nsoilay)      ! wilting soil moisture value (fraction of pore space)
!     rhosoi(npoi,nsoilay)     ! soil density (without pores, not bulk) (kg m-3)
!     poros(npoi,nsoilay)      ! porosity (mass of h2o per unit vol at sat / rhow)
!     porosflo(npoi,nsoilay)   ! porosity after reduction by ice content
!     sand(npoi,nsoilay),      ! percent sand of soil
!     clay(npoi,nsoilay),      ! percent clay of soil
!     stressl(npoi,nsoilay)    ! soil moisture stress factor for the  lower canopy (dimensionless)
!     stressu(npoi,nsoilay)    ! soil moisture stress factor for the upper canopy (dimensionless)
!     upsoiu(npoi,nsoilay)     ! soil water uptake from transpiration (kg_h2o m-2 s-1)
!     upsoil(npoi,nsoilay)     ! soil water uptake from transpiration (kg_h2o m-2 s-1)

      real*8, dimension(:,:), allocatable :: cpwf

!     cpwf(npoi, nsoilay)      ! Kaiyuan Li for Green-Ampt, capillary presssure at wetting front (m)
!     swater(npoi, nsoilay)    ! Kaiyuan Li for Green-Ampt, to save the value of wsoi before infiltration
!     sice(npoi, nsoilay)      ! Kaiyuan Li for Green-Ampt, to save the value of wisoi before infiltration

       real*8, dimension(:), allocatable :: stre_tl, stre_tu, fwpudtot

!     stre_tl(npoi)            ! Kaiyuabn Li for water uptake 
!     stre_tu(npoi)            ! Kaiyuan Li for water uptake  
!     fwpudtot(npoi)           ! Kaiyuan Li for Green-Ampt, accumulated infiltration during a single rain event (mm)

!     tnuptake(npoi,nsoilay)   ! total potential nitrogen uptake for current timestep (kg_n/timestep)
!     anuptake(npoi,nsoilay)   ! actual nitrogen uptake for current timestep (kg_n/timestep) 
!     smsoil(npoi,-nsoilay:nsoilay)    ! current timestep solute in soil (kg_solute/m2)
!     smsoln(npoi,-nsoilay:nsoilay)    ! current timestep solute in solution (kg_solute/m2)
!     fout(npoi,-nsoilay:nsoilay)      ! current timestep total nitrogen leaching (flux)  (kg_solute/m2)
!     nout(npoi,-nsoilay:nsoilay)      ! current timestep nitrate leaching flux  (kg_solute/m2)
!     csoln(npoi,-nsoilay:nsoilay)     ! current timestep solute concentration in solution (mg/liter)
!     domtext(npoi,nsoilay),
!     fracsand(npoi,nsoilay),
!     fracsilt(npoi,nsoilay),
!     fracclay(npoi,nsoilay),
!     swater(npoi, nsoilay),   ! Green-Ampt, to save the value of wsoi before infiltration
!     sice(npoi, nsoilay)      ! Green-Ampt, to save the value of wisoi before infiltration
! 
      real*8, dimension(:,:), allocatable :: hflo
!
!     hflo(npoi,nsoilay+1)     ! downward heat transport through soil layers (W m-2)
!
      integer, dimension(:,:), allocatable :: ibex
!
!     ibex(npoi,nsoilay)       ! nint(bex), used for cpu speed
!
      real*8, dimension(:,:), allocatable :: qglif
!
!     qglif(npoi,4)            ! 1: fraction of soil evap (fvapg) from soil liquid
!                              ! 2: fraction of soil evap (fvapg) from soil ice
!                              ! 3: fraction of soil evap (fvapg) from puddle liquid
!                              ! 4: fraction of soil evap (fvapg) from puddle ice

      integer nslaym           ! number of soil layers to 1 m depth

      real*8, dimension(:,:), allocatable :: wflo
!     wflo(npoi,nsoilay+1),    !downward h20 flow across boundaries  

      real*8, dimension(:), allocatable :: fwtop, fwpud, deposn
!     fwtop(npoi),            ! h20 flux into top layer from evap/condens 
!     fwpud(npoi),            ! h20 flux into top layer from puddle (kg m-2/s)
!     deposn(npoi)            ! daily nitrogen deposition from atmosphere (kg n m-2 day-1)
!
      real*8 co, cf, fnitrate, nloss, cndepth
!      co             ! initial inorganic nitrogen in soil storage pools 
!      cf             ! buffering constant for inorganic N in soil solution
!      fnitrate       ! fraction of leached inorganic nitrogen that is nitrate - other portion is ammonium
!      nloss          ! fraction of immobile inorganic nitrogen that is lost annually to volatilization/denitrification 



!
end module inland_comsoi
