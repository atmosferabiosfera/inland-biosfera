#include "inland_config.h"
module inland_comsatp

      implicit none
      public
      save

! ------
! comsatp
! ------
!
! ---------------------------------------------------------------------
! parameters associated with statement functions
! ---------------------------------------------------------------------
!
! polynomials for svp(t), d(svp)/dt over water and ice are from
! lowe(1977),jam,16,101-103.
      real*8 asat0, asat1, asat2, asat3, asat4, asat5, asat6

      parameter (asat0 =  6.1078000, &
                asat1 =  4.4365185e-1, &
                asat2 =  1.4289458e-2, &
                asat3 =  2.6506485e-4, &
                asat4 =  3.0312404e-6, &
                asat5 =  2.0340809e-8, &
                asat6 =  6.1368209e-11 )

      real*8 bsat0, bsat1, bsat2, bsat3, bsat4, bsat5, bsat6

      parameter (bsat0 =  6.1091780, &
                bsat1 =  5.0346990e-1, &
                bsat2 =  1.8860134e-2, &
                bsat3 =  4.1762237e-4, &
                bsat4 =  5.8247203e-6, &
                bsat5 =  4.8388032e-8, &
                bsat6 =  1.8388269e-10 )

      real*8 csat0, csat1, csat2, csat3, csat4, csat5, csat6

      parameter (csat0 =  4.4381000e-1, &
                csat1 =  2.8570026e-2, &
                csat2 =  7.9380540e-4, &
                csat3 =  1.2152151e-5, &
                csat4 =  1.0365614e-7, &
                csat5 =  3.5324218e-10, &
                csat6 = -7.0902448e-13 )

      real*8 dsat0, dsat1, dsat2, dsat3, dsat4, dsat5, dsat6

      parameter (dsat0 =  5.0303052e-1, &
                 dsat1 =  3.7732550e-2, &
                 dsat2 =  1.2679954e-3, &
                 dsat3 =  2.4775631e-5, &
                 dsat4 =  3.0056931e-7, &
                 dsat5 =  2.1585425e-9, &
                 dsat6 =  7.1310977e-12 )
!
end module inland_comsatp
