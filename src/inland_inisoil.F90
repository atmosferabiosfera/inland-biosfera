#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine inisoil (irestart)
! ---------------------------------------------------------------------
! does initialization for soil database
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsoi
      use inland_comsno
      use inland_comtex
      use inland_comsum
      use inland_comcrop, only:isimagro, totirrig
      use inland_comnitr
      implicit none
!------------------------------Arguments--------------------------------
! Input arguments
      integer irestart     ! 0: not restart of INLAND, 1 restart of INLAND

! local variables
      integer i,       &
              k,       &
              l,       &
              msand,   &       ! % of sand in grid point
              mclay,   &       ! % of clay in grid point
              lmin,    &       ! closest textural class from texture of grid point
              textcls, &       ! textural class assignment (1..11)
              domtexti

      real*8 fsand,  &       ! fraction of sand in grid point
             fsilt,  &       ! fraction of silt in grid point
             fclay,  &       ! fraction of clay in grid point
             forganic        ! fraction of organic matter in soil

      real*8 xdat(ndat), &
             ydat(ndat), &
             zdat(ndat), &
             dsq(ndat)

     if(isimagro .gt. 0) then
! --------------------------------------------------------------------------
! added for Green-Ampt infiltration model
! capillary pressure at the wetting front for Green-Ampt infiltration model
! units are meters of head (from Rawls et al. 1983, ASCE 109(1), Table 2)
!
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! added for Green-Ampt
      do l = lbeg, lend
         fwpudtot(l) = 0.0
        
        if (irestart .eq. 0) then
         wpud(l) = 0.0
        endif
      end do
!
     endif ! check for crop existence
! -----------------------------------------------------------------------
! set sand/silt/clay vectors (xdat,ydat,zdat) for 11 data points
      do 100 l = 1, ndat
         xdat(l) = texdat(1,l)
         ydat(l) = texdat(2,l)
         zdat(l) = texdat(3,l)
100   continue

! initialization and normalization constant for puddle model (kg m-2)
      !wpudmax = 4.5
      if (irestart .eq. 0) then
         wpud(:)=0.0
         wipud(:)=0.0
         fwpudtot(:)=0.0
      endif

! set prescribed soil layer thicknesses
!      hsoi(1) = 0.10     hsoi(1)  = 0.05
!      hsoi(2) = 0.15     hsoi(2)  = 0.05
!      hsoi(3) = 0.25     hsoi(3)  = 0.1
!      hsoi(4) = 0.50     hsoi(4)  = 0.1
!      hsoi(5) = 1.00     hsoi(5)  = 0.1
!      hsoi(6) = 2.00     hsoi(6)  = 0.2
!                         hsoi(7)  = 0.2
!                         hsoi(8)  = 0.2
!                         hsoi(9)  = 0.5
!                         hsoi(10) = 0.5
!                         hsoi(11) = 0.5

! set physical parameters of soil
      z0soi(:)=0.005

! initialize soil water and soil temperature fields
      if (irestart .eq. 0) then

         wisoi(:,:)=0.00

! wsoi is now read from params/soil
#ifndef SINGLE_POINT_MODEL
         tsoi(:,:)=278.13
     if(isimagro .gt. 0) then
         wsoi(:,:)=0.50
         adwsoilay(:,:)=0.50
         adwisoilay(:,:)=0.00
!         adtsoilay(:,:)=278.13
! initialize total irrigation
         totirrig(:)=0.00
! initialize plant available nitrogen pool (kg/m2)  
         aplantn(:)=0.0080
! initialize soil  nitrogen variables
         totnuptake(:,:)=0.00
         stressn(:,:)=1.0

     endif ! check for crop existence

         tg(:)=278.13
         ti(:)=273.13
#else /* SINGLE_POINT_MODEL */
! FIXME: This change is not documented anywhere. - fzm
! tsoi for single point is read from single_point_parameters
         tg(:)=297.66
         ti(:)=273.16
#endif /* SINGLE_POINT_MODEL */
      endif

! set soil surface parameters for the global domain
     if(isimagro .eq. 0) then

      do 200 i = lbeg, lend
!
! SEE BELOW OTHER PARAMETERS
!
! Convert input sand and clay percents to fractions
         msand = nint(sand(i,1))
         mclay = nint(clay(i,1)) 
         fsand = 0.01 * msand
         fclay = 0.01 * mclay
         fsilt = 0.01 * (100 - msand - mclay)

! soil surface albedo:
!
! from bats table 3.ii assuming albedo depends on texture
#ifndef SINGLE_POINT_MODEL
         albsav(i) = fsand * 0.120 + fsilt * 0.085 + fclay * 0.050
         albsan(i) = 2.0 * albsav(i)
#else /* SINGLE_POINT_MODEL */
! FIXME: There's no explanation on why the code turned into this for 0D version
!       - fzm
         albsav(i) = 0.10
         albsan(i) = 0.40
#endif /* SINGLE_POINT_MODEL */
200 continue
    endif ! check for crop existence


! create soil properties look-up table
!
! set soil parameters at each layer for the global domain
! soita.nc file is for layers only to 4 m; currently this
! is for a total of six layers.  For the remaining  
! layers below that, set texture to be equal to that of the
! last layer (layer 6)
! analysis of the current WISE-IGBP soil textural dataset
! reveals very little information below 4 m.

! gabriel abrahao: in order to make natural and agro more compatible, 
! the default input files now have 11 layers, ending at 240cm. 
! The information for the 240cm layer is the same as the 400cm layer
! in the old default. Also, both parameter files have 12 layers
! Its just a matter of repeating the 11th one for layers below it,
! and the natural vegetation will be using essentially the same information
! as before but solving the equations with more layers. 
      do 310 i = lbeg, lend 
         do 300 k = 1, nsoilay

! Convert input sand and clay percents to fractions
            !Fontes: infilensoilayer is equal 0 on single_point
            if (k.le.infilensoilayer.or.infilensoilayer.eq.0) then
               msand = nint(sand(i,k))
               mclay = nint(clay(i,k)) 
            else
               msand = nint(sand(i,infilensoilayer)) 
               mclay = nint(clay(i,infilensoilayer)) 
          endif
      if(isimagro .eq. 0) then
          fsand    = 0.01 * msand
          fclay    = 0.01 * mclay
          fsilt    = 0.01 * (100 - msand - mclay)
     endif


! for now, we assume that all soils have a 1% organic content -- 
! this is just a place holder until we couple the soil carbon
! dynamics to the soil physical properties
            forganic = 0.010

! density of soil material (without pores, not bulk) (kg m-3)
! from Campbell and Norman, 1998
            rhosoi(i,k) = 2650.0 * (1.0 - forganic) + 1300.0 * forganic 

! specific heat of soil material (j kg-1 k-1):
! from Campbell and Norman, 1998
            csoi(i,k) =  870.0 * (1.0 - forganic) + 1920.0 * forganic 

! C. Kucharik
! match textural fractions with soil textural class 
! calls two functions to match sand and clay fractions
! with proper soil textural class based on the usda
! classification system
!
! Now use input dominant soil texture class from soil_text.nc
! Must convert from CONUS texture class values to those of IBIS
! no silt in IBIS so set CONUS silt equal to silt loam
!
          lmin = textcls (msand,mclay)  !class from the global file.



          fracsand(i,k) = texdat(1,lmin)
          fracsilt(i,k) = texdat(2,lmin)
          fracclay(i,k) = texdat(3,lmin)

! Kaiyuan Li for Green-Ampt
            cpwf(i, k) = cpwfdat(lmin)
!           swater(i, k) = 0.01   ! save wsoi before infiltration
            swater(i, k) = 0.00000001 
            sice(i, k) = 0        ! save wisoi before infiltration

! porosity (fraction):
            poros(i,k) = porosdat(lmin)

! field capacity (defined relative to the porosity):
            sfield(i,k) = 1.0 / poros(i,k) * sfielddat(lmin)

! wilting point (defined relative to the porosity):
            swilt(i,k)  = 1.0 / poros(i,k) * swiltdat(lmin)

! "b" exponent for the Campbell moisture-release equation:
            bex(i,k) = bexdat(lmin)

! nearest integer of "b" exponent (for computational efficiency):
            ibex(i,k) = nint(bex(i,k))

! saturated matric (air entry) potential (m-h2o):
            suction(i,k) = suctiondat(lmin)

! saturated hydraulic conductivity (m s-1):
          hydraul(i,k) = hydrauldat(lmin)

! -----------------------------------------------------------------------
! added for Green-Ampt
      if(isimagro .gt. 0)then
          swater(i, k) = 0.000001   ! useless value
          sice(i, k) = 0        ! useless value
      endif
! -----------------------------------------------------------------------

 300  continue

! surface parameters
      if(isimagro .gt. 0)then
        albsav(i) = fracsand(i,1) * 0.120 + &
                    fracsilt(i,1) * 0.085 + &
                    fracclay(i,1) * 0.050
        albsav(i) = 0.07
        albsan(i) = 2.0 * albsav(i)
      endif

310  continue
!
      return
end subroutine inisoil
