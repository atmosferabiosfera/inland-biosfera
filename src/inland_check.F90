#include "inland_config.h"
#ifndef SINGLE_POINT_MODEL
subroutine check(irestart, soilcspin, nrun, snorth, ssouth, swest, seast)
#else /* SINGLE_POINT_MODEL */
subroutine check(soilcspin, nrun, isoilforc)
#endif /* SINGLE_POINT_MODEL */

! ---------------------------------------------------------------------
      use inland_control, only: iyearout, imonthout, idailyout, isimveg, &
                                isimfire, isimco2, isinfilt,isimrwu, isimland, &
                                itauw, ivmax, isla, ica
      implicit none

      integer soilcspin,   & 
              irestart,    &
	      nrun
      real*8 snorth,       &
             ssouth,       &
             swest,        &
             seast
#ifdef SINGLE_POINT_MODEL
      integer isoilforc
#endif /* SINGLE_POINT_MODEL */

#ifndef SINGLE_POINT_MODEL
! verify that domain values were read and are correct
      if ( snorth .eq. 999.0 .or. ssouth .eq. 999.0 .or. swest .eq. 999.0 .or. seast .eq. 999.0 ) then
         write (*,*) 'ERROR: domain dimensions not set'
         stop 1
      end if
      if ( snorth .gt. 90.0 .or. snorth .lt. -90.0 .or. ssouth .gt. 90.0 .or. ssouth .lt. -90.0 .or. &
           swest .gt. 180.0 .or. swest .lt. -180.0 .or. seast .gt. 180.0 .or. seast .lt. -180.0 ) then
         write (*,*) 'ERROR: domain dimensions exceed bounds (+90N -90S -180W +180E)'
         write (*,*) snorth, ssouth, swest, seast
         stop 1
      end if
      if ( snorth .lt. ssouth .or. swest .gt. seast ) then
         write (*,*) 'ERROR: domain dimensions must be N>S / W<E'
         write (*,*) snorth, ssouth, swest, seast
         stop 1
      end if

      if ( itauw .lt. 0 .or. ivmax .lt. 0 .or. isla .lt. 0 .or. ica .lt. 0 .or. &
           itauw .gt. 1 .or. ivmax .gt. 1 .or. isla .gt. 1 .or. ica .gt. 1) then
         write(*,*) 'ERROR: itauw, ivmax, isal and ica must be between 0 and 1'
         write(*,*) 'itauw: ', itauw, 'ivmax: ',ivmax,'isla:' ,isla,'ica:' ,ica
         stop 1          	 
      end if

! Check output parameters
      if ( iyearout .lt. 0 .or. imonthout .lt. 0 .or. idailyout .lt. 0 .or.  &
           iyearout .gt. 1 .or. imonthout .gt. 1 .or. idailyout .gt. 1) then
         write(*,*) 'ERROR: Values of output parameters less 0 or exceed 1'
         write(*,*) 'iyearout: ', iyearout, 'imonthlyout: ',imonthout,'idailyout:' ,idailyout
         stop 1          	 
      end if

! Check land use parameters
      if ( isimland .lt. 0 .or. isimland .gt. 1) then
         write(*,*) 'ERROR: isimland must be between 0 and 1'
         write(*,*) 'isimland: ' ,isimland
         stop 1          	 
      end if


#endif /* SINGLE_POINT_MODEL */
      
! Check simulation parameters
      if ( soilcspin .lt. 0 .or. isimco2 .lt. 0 .or. soilcspin .gt. 1 .or. &
           isimco2 .gt. 1) then
         write(*,*) 'ERROR: Values of simulation parameters must be between 0 and 1'
         write(*,*) 'isimco2: ', isimco2, 'soilcspin: ', soilcspin 
         stop 1          	 
      end if
      
      if ( isimfire .lt. 0 .or. isimfire .gt. 3) then
         write(*,*) 'ERROR: isimfire must be between 0 and 3'
         write(*,*) 'isimfire: ' ,isimfire
         stop 1          	 
      end if
      
      if ( isimveg .lt. 0 .or. isimveg .gt. 2) then
         write(*,*) 'ERROR: isimveg must be between 0 and 2'
         write(*,*) 'isimveg: ' ,isimveg
         stop 1          	 
      end if

      if ( isinfilt .lt. 0 .or. isinfilt .gt. 1) then
         write(*,*) 'ERROR: isinfilt must be between 0 and 1'
         write(*,*) 'isinfilt: ' ,isinfilt
         stop 1          	 
      end if

      if ( isimrwu .lt. 0 .or. isimrwu .gt. 1) then
         write(*,*) 'ERROR: isimrwu must be between 0 and 1'
         write(*,*) 'isimrwu: ' ,isimrwu
         stop 1          	 
      end if
 
#ifdef SINGLE_POINT_MODEL      
      if ( isoilforc .lt. 0 .or. isoilforc .gt. 1) then
         write(*,*) 'ERROR: isoilforc must be between 0 and 1'
         write(*,*) 'isoilforc: ' ,isoilforc
         stop 1          	 
      end if
#endif /* SINGLE_POINT_MODEL */      

! Check temporal parameters

#ifndef SINGLE_POINT_MODEL
      if ( irestart .lt. 0 .or. irestart .gt. 1) then
         write(*,*) 'ERROR: irestart must be between 0 and 1'
         write(*,*) 'irestart: ', irestart
         stop 1          	 
      end if
#endif /* SINGLE_POINT_MODEL */

      if ( nrun .lt. 0) then
         write(*,*) 'ERROR: Value of nrun less 0'
         write(*,*) 'nrun :', nrun
         stop 1          	 
      end if

      return
end subroutine check
