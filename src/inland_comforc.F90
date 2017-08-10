#include "inland_config.h"
! inland_comforc.h - Ported from Hewlley Imbuzeiro's 0D INLAND version

! comfor.h last update 11/11/06 S.V. Cuadra
!
! This module holds the parameters and for the forcing subroutines from 0D model

! -------
! comforc
! -------

#ifndef SINGLE_POINT_MODEL
#error "This module should ONLY be compiled for 0D INLAND model option."
#endif /* SINGLE_POINT_MODEL */

module inland_comforc
      implicit none
      public
      save

      integer dimforc                        ! dimension (number of lines) of input data

!            xlati(npoi)                     ! latitude of the considered point
      real*8, dimension(:), allocatable :: xlati

!            xta(npoi,dimforc),            & ! air temperature at 64m (K)
!            xprec(npoi,dimforc),          & ! total precipitation at 64m (mm)
!            xua(npoi,dimforc),            & ! wind speed at 64m (m.s-1)
!            xcld(npoi,dimforc),           & ! cloudiness (%)
!            xqa(npoi,dimforc),            & ! relative humidity at 64m 
!                                             (kg_h2o/kg_air)
!            xsin(npoi,dimforc),           & ! incoming solar radiation (W/m-2)
!            xlin(npoi,dimforc),           & ! incoming longwave radiatn.(W/m-2)
      real*8, dimension(:,:), allocatable :: xta, xprec, xua, xcld, xqa, xsin, &
                                             xlin

!            xsmoi(npoi,nsoilay,dimforc),  & ! soil moisture (fraction)
!            xstemp(npoi,nsoilay,dimforc), & ! soil temperature (K)
      real*8, dimension(:,:,:), allocatable :: xsmoi, xstemp
end module inland_comforc
