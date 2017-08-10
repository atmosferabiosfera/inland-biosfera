#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine vdiag
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control, only: iyear, iyear0
      use inland_comveg
      use inland_comsum
      use inland_subgrid

      implicit none
!------------------------------Variables--------------------------------
! local variables
      integer i, k, j, ilpt

      real*8 vtarea(nvegtype),   & ! total area of the vegetation type (m**2)
             vnee(nvegtype),     & ! vegetation type average nee (kg-c/m**2/yr)
             vnpp(nvegtype),     & ! vegetation type average npp (kg-c/m**2/yr)
             vgpp(nvegtype),     & ! vegetation type average gpp (kg-c/m**2/yr)
             vleach(nvegtype),   & ! vegetation type average leaching (kg-c/m**2/yr)
             vco2mic(nvegtype),  & ! vegetation type average micr. respiration (kg-c/m**2/yr)
             vcdist(nvegtype),   & ! vegetation type average disturbance (kg-c/m**2/yr)
             vtotfall(nvegtype), & ! vegetation type average litterfall (kg-c/m**2/yr)
             vbiomass(nvegtype), & ! vegetation type average biomass (kg-c/m**2)
             vlai(nvegtype),     & ! vegetation type average lai (m**2/m**2)
             vsoic(nvegtype),    & ! vegetation type average soil carbon (kg-c/m**2)
             vrunoff(nvegtype)     ! vegetation type average runoff (mm/yr)

      real*8 scaling          ! scaling to apply to each point = garea(i) * tilefrac(i)


! initialize variables
      do 100 k = 1, nvegtype
         vnee(k)     = 0.0
         vtarea(k)   = 0.0
         vnpp(k)     = 0.0
         vgpp(k)     = 0.0
         vbiomass(k) = 0.0
         vlai(k)     = 0.0
         vsoic(k)    = 0.0
         vrunoff(k)  = 0.0
         vleach(k)   = 0.0
         vco2mic(k)  = 0.0
         vcdist(k)   = 0.0
         vtotfall(k) = 0.0
100   continue

! sum ecosystem properties over each vegetation type
      do j = 1, npoi1
         do ilpt = 1, mlpt

            i = subgrid_get_index(j,ilpt)
            ! for now garea is the area of the entire cell, if this changes must remove tilefrac here
            scaling = scaling * tilefrac(i)

         k = int (max (dble(1.0), min (dble(nvegtype), vegtype0(i))))
         vtarea(k) = vtarea(k)    +scaling
         vnee(k)     = vnee(k)    +scaling*ayneetot(i)
         vnpp(k)     = vnpp(k)    +scaling*aynpptot(i)
         vgpp(k)     = vgpp(k)    +scaling*aygpptot(i)
         vbiomass(k) = vbiomass(k)+scaling*totbiou(i)+scaling*totbiol(i)
         vlai(k)     = vlai(k)    +scaling*totlaiu(i)+scaling*totlail(i)
         vsoic(k)    = vsoic(k)   +scaling*aycsoi(i) 
         vrunoff(k)  = vrunoff(k) +scaling*aytrunoff(i)
         vleach(k)   = vleach(k)  +scaling*yrleach(i)
         vco2mic(k)  = vco2mic(k) +scaling*ayco2mic(i)
         vcdist(k)   = vcdist(k)  +scaling*cdisturb(i)
         vtotfall(k) = vtotfall(k)+scaling*(falll(i)+fallr(i)+fallw(i))
         end do ! mlpt
      end do ! npoi1

!
! calculate area averages
      do 300 k = 1, nvegtype
         vnee(k)     = vnee(k)     / max (dble(1.0), vtarea(k))
         vnpp(k)     = vnpp(k)     / max (dble(1.0), vtarea(k))
         vgpp(k)     = vgpp(k)     / max (dble(1.0), vtarea(k))
         vbiomass(k) = vbiomass(k) / max (dble(1.0), vtarea(k))
         vlai(k)     = vlai(k)     / max (dble(1.0), vtarea(k))
         vsoic(k)    = vsoic(k)    / max (dble(1.0), vtarea(k))
         vrunoff(k)  = vrunoff(k)  / max (dble(1.0), vtarea(k))
         vleach(k)   = vleach(k)   / max (dble(1.0), vtarea(k))
         vco2mic(k)  = vco2mic(k)  / max (dble(1.0), vtarea(k))
         vcdist(k)   = vcdist(k)   / max (dble(1.0), vtarea(k))
         vtotfall(k) = vtotfall(k) / max (dble(1.0), vtarea(k))
300   continue

! write some diagnostic output to history file (file opened in inland_main_offline.F90)
! TODO: wasnt open file=202
      if (myid .eq. 0) then
         do 400 k = 1, nvegtype
         write (202,9000) iyear, k, vtarea(k)/1.0e+06, vnee(k), vnpp(k), &
                          vgpp(k), vleach(k), vco2mic(k), vcdist(k),     &
                          vtotfall(k), vbiomass(k), vlai(k), vsoic(k),   &
                          vrunoff(k)
400      continue
!        call flush (202)
      endif
!
9000  format (1x,i4,5x,i2,5x,1e10.3,11f10.3)
!
! return to main program
!
      return
end subroutine vdiag
