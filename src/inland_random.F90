#include "inland_config.h"

real function ran2(idum,idum2,iv,iy)
      use inland_control, only: env_ran2val
      implicit none
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real am, eps, rnmx
      parameter (im1=2147483563,   &
                 im2=2147483399,   &
                 am=1./im1,        &
                 imm1=im1-1,       &
                 ia1=40014,        &
                 ia2=40692,        &
                 iq1=53668,        &
                 iq2=52774,        &
                 ir1=12211,        &
                 ir2=3791,         &
                 ntab=32,          &
                 ndiv=1+imm1/ntab, &
                 eps=1.0e-7,       &
                 rnmx=1.-eps)

      integer idum2,j,k,iv(ntab),iy

      if ( env_ran2val .gt. 0. ) then
         ran2=env_ran2val
         return
      end if

      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 10 j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum=idum+im1
            if (j.le.ntab) iv(j)=idum
10       continue
         iy=iv(1)
      endif

      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1

      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2

      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1

      ran2=min(am*iy,rnmx)
      return
end function ran2

!*****************************************************************************80
!
!! R4_UNI returns a uniformly distributed real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Output, real ( kind = 4 ) R4_UNI, a uniformly distributed random value in
!    the range [0,1].
!
! taken from
! http://people.sc.fsu.edu/~jburkardt/f_src/ziggurat_openmp/ziggurat_openmp.f90
!

function r4_uni ( jsr )

  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_input
  real ( kind = 4 ) r4_uni

  jsr_input = jsr

  jsr = ieor ( jsr, ishft ( jsr,   13 ) )
  jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
  jsr = ieor ( jsr, ishft ( jsr,    5 ) )

! r4_uni = 0.5E+00 + 0.2328306E-09 * real ( jsr_input + jsr, kind = 4 )

  r4_uni = 0.5E+00 + real ( jsr_input + jsr, kind = 4 ) &
    / real ( 65536, kind = 4 ) / real ( 65536, kind = 4 )

  return
end


! this function used by inland

function ran_uni ( seed )

  use inland_control, only: env_ran2val
  implicit none

  integer seed
  real ( kind = 4 ) ran_uni
  real r4_uni, ran2

  ! depending on INLAND_RANDOMVAL value, use different algorithm
  ! 0 (default) : r4_uni, recommended for thread-safety, with a seed for each grid point
  ! > 0         : use fixed value
  ! < 0         : call random_number, not thread-safe
  if ( env_ran2val .eq. 0. ) then
     ran_uni = r4_uni( seed )
  else if ( env_ran2val .gt. 0. ) then
     ran_uni = env_ran2val
  else ! if ( env_ran2val .lt. 0. )
     call random_number(ran_uni)
  end if

  return
end
