#include "inland_config.h"
! ------------------------------------------------------------------------
      real*8 function twet3(tak, qq1, pp1)
! ------------------------------------------------------------------------
!
! twet3.f last update 8/30/2000 C Molling
!
! This function calculates the wet bulb temperature given
! air temp, specific humidity, and air pressure.  It needs the function esat
! in order to work (in comsat.h).  The function is an approximation to
! the actual wet bulb temperature relationship.  It agrees well with the
! formula in the Smithsonian Met. Tables for moderate humidities, but differs
! by as much as 1 K in extremely dry or moist environments.
!
! INPUT
!     tak - air temp in K
!     qq1 - specific humidity in kg/kg			!Castanho Kai Li for g77 p, q -> pp1, qq1
!     pp1 - air pressure in Pa (Pa = 100 * mb)		!Castanho Kai Li for g77 p, q -> pp1, qq1
!
! OUTPUT
!     twet3 - wet bulb temp in K, accuracy?
!
      use inland_parameters
      use inland_comsatp

      implicit none

!------------------------------Arguments--------------------------------
      real*8 tak   ! air temp in K
      real*8 qq1     ! specific humidity in kg/kg
      real*8 pp1     ! air pressure in Pa (Pa = 100 * mb)

!------------------------------Variables--------------------------------
! local variables
      integer i

      real*8 ta, twk, twold, diff

#define TWET3_COMSAT
#include "inland_comsat.h"

! temperatures in twet3 equation must be in C
! pressure in qsat function must be in Pa
! temperatures in esat,hvapf functions must be in K
!
!     Air temp in C
!     -------------
      ta = tak - 273.16
!
!     First guess for wet bulb temp in C, K
!     -------------------------------------
      twet3 = ta * qq1 / qsat(esat(tak),pp1)
      twk = twet3 + 273.16
!
!     Iterate to converge
!     -------------------
      do 100 i = 1, 20
         twold = twk - 273.16
         twet3 = ta - (hvapf(twk,tak)/cair) * ( qsat( esat(twk),pp1 )-qq1 )
         diff = twet3 - twold
!
! below, the 0.2 is the relaxation parameter that works up to 40C (at least)
!
         twk = twold + 0.2 * diff + 273.16
         if (abs(twk-273.16-twold) .lt. 0.02) goto 999
 100  continue
!
      if (myid .eq. 0) then
         print *, 'Warning, twet3 failed to converge after 20 iterations!'
         print *, 'twet3, twetold: ', twk, twold+273.16
         print *, 'twetbulb is being set to the air temperature'
      endif
!
      twet3 = tak
!
!     Return wet bulb temperature in K
!     --------------------------------
 999  twet3 = twk
!
      end function
