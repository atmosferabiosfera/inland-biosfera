#include "inland_config.h"

!   This subroutine Allocates variables as they are defined by rd_param.

subroutine inland_inneralloc

    use inland_comsoi, only: hsoi
    use inland_parameters, only: npft,nband,nsoilay,nlon,nlat,nsnolay
    use inland_comwork, only: ndim2,ndim3,ndim4,work,cdummy,cdummyint

    ! Allocate variables in comsoi
    allocate(hsoi(nsoilay+1))! FIXME: really needed to allocate nsoilay+1?
    hsoi(:) = 0.

    ! Get comwork's ndim4 size
!    The base values are now in inland_parameters.F90 module; 
!    in the future they will be read from inland.infile
    ndim2 = nlon*nlat
    ndim4 = max(nlon,nlat,nband,nsoilay,nsnolay,npft)
    ndim3=nlon*nlat*max(nband,nsoilay,nsnolay,npft)
    if ( mlpt .gt. 1 ) ndim3 = ndim3*(mlpt+1)

    allocate(work(ndim2))
    allocate(cdummy(ndim3),cdummyint(ndim3))
    work(:) = 0.
    cdummy(:) = 0.
    cdummyint(:) = 0


end subroutine inland_inneralloc
