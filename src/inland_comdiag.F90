#include "inland_config.h"
! inland_comdiag.F90 - last update 2011-06-28 - FZM
! comdiag.h last update was 2/12/01 C. Molling

! This holds the parameters and common blocks for diagnostics
module inland_comdiag
      implicit none
      public
      save

! TODO: this diagnostic stuff may be dynamically allocated and read from
! files. - fzm
      integer nvars, nfiles
      parameter (nvars=128, nfiles=10)   ! nvars = # diagnostic variables

      integer ldiag(nvars,nfiles)      ! chosen diagnostic variables
!      common /diag1/ ldiag

      integer diagstart(nfiles),    & ! year diagnostic output begins
              diagend(nfiles),      & ! year diagnostic output stops
              ndiagpt(nfiles),      & ! point in an npoi array
              nfreq(nfiles)           ! frequency of diagnostic output

!      common /diag2/diagstart,diagend,ndiagpt,nfreq
end module inland_comdiag
