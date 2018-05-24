#include "inland_config.h"
module inland_comatm

      implicit none
      public
      save
!
! ------
! comatm
! ------
!
      real*8, dimension(:), allocatable :: coszen, fira, solsoi
!
!     coszen(npoi)       ! cosine of solar zenith angle
!     fira(npoi)         ! incoming ir flux (W m-2)
! \/\/\/\/\/\/ ALLOC!!!!
!     solsoi(npoi),      ! downward solar flux at soil (W m-2)
!
      real*8, dimension(:,:), allocatable :: solad, solai, asurd, asuri
!
!     solad(npoi,nband)  ! direct downward solar flux (W m-2)
!     solai(npoi,nband)  ! diffuse downward solar flux (W m-2)
!     asurd(npoi,nband)  ! direct albedo of surface system
!     asuri(npoi,nband)  ! diffuse albedo of surface system 


!
      real*8, dimension(:), allocatable :: ua, ux, uy, ta, qa, raina, snowa
!
!     ua(npoi)           ! wind speed (m s-1)
!     ux(npoi)           ! northward wind (m s-1)
!     uy(npoi)           ! eastward wind (m s-1)
!     ta(npoi)           ! air temperature (C)
!     qa(npoi)           ! specific humidity (kg_h2o/kg_air)
!     raina(npoi)        ! rainfall rate (mm/day)
!     snowa(npoi)        ! snowfall rate (mm/day)
!
      real*8, dimension(:), allocatable :: psurf, tmax, tmin, qd, ud
!
!     psurf(npoi)        ! surface pressure (Pa)
!     tmax(npoi)         ! maximum daily temperature (K)
!     tmin(npoi)         ! maximum daily temperature (K)
!     qd(npoi)           ! daily average specific humidity (kg_h2o/kg_air)
!     ud(npoi)           ! daily average wind speed (m/sec)
!
      real*8 co2conc, o2conc
!
!     co2conc            ! co2 concentration (mol/mol)
!     o2conc             ! o2 concentration (mol/mol)

!--------------------------------------------------------------------------------
! Castanho Kaiyuan Li for Green-Ampt

      real*8 t_sec, t_startp

!     t_sec,     ! time since midnight in second
!     t_startp   ! time since midnight in second when precipitation starts 
! -------------------------------------------------------------------------------

! Weather generator-specific variables
!---------------------------------------------
! FIXME: if these variables are used only by inland_daily subroutine, why the
!       need to make part of a module?
!       Common variables on daily and diurnal subroutines: precip and cloud.
!       ALL other variables are used only on inland_daily subroutine.
! precip(npoi): daily precipitation (mm/day)
! precipdaysum(npoi): cumulaive precipitation in mm for the month
! cloud(npoi): cloud fraction
      real*8,  dimension(:), allocatable :: precip,precipdaysum,cloud

! precipday(npoi,31): precipitation in mm for a give dy of the month
! xstore(npoi,3): weather generator 'memory' matrix FIXME: better variable name
      real*8,  dimension(:,:), allocatable :: precipday,xstore

! iwet(npoi): wet day / dry day flag
! iwetdaysum(npoi): counter for wet days in the current month
      integer, dimension(:), allocatable :: iwet,iwetdaysum

! iwetday(npoi,31): days that were wet along the month
      integer, dimension(:,:), allocatable :: iwetday

! FIXME: if these variables are used only by inland_daily subroutine, why the
!       need to make part of a module?
      real*8, dimension(:), allocatable :: rh
end module inland_comatm
