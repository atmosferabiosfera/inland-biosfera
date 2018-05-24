#include "inland_config.h"
!---------------------------------------------------------------------
      subroutine inidiag(idiag)
!---------------------------------------------------------------------
! CALLED BY main
! CALLS diaginit
! This subroutine reads in diagnostics info from diag.infile and calls
! the subroutine diaginit.
!
! INPUT
      integer idiag          ! number of diagnostic files
! COMMON
      use inland_comdiag.F90
! INTERNAL
      integer lun, i, ivar
!
      real*8  diaglat(nfiles),&      ! latitude of diagnostic point
              diaglon(nfiles)        ! longitude of diagnostic point
!
        lun = 12
        open (lun, status = 'old', file = 'diag.infile')
!
        do 100 i = 1, 26
          read (lun,*)
 100    continue
!
        do 110 i = 1, nfiles
          read (lun,*) diaglat(i)
          read (lun,*) diaglon(i)
          read (lun,*) diagstart(i)
          read (lun,*) diagend(i)
          read (lun,*) nfreq(i)
          read (lun,*)
 110    continue
!
        do 120 i = 1, 8
          read (lun,*)
 120    continue
!
        do 130 ivar = 1, nvars
          read (lun,*) (ldiag(ivar,ifile),ifile=1,nfiles)
 130    continue
!
        close (lun)
!
         call diaginit (idiag, diaglat, diaglon)
!
         return
end subroutine inidiag
