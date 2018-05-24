#include "inland_config.h"
! map_offline: map coupled model variables used in the model given
! the input data.

subroutine map_offline(imonthp, rdlsf, loopi, kpti, kptj, ipointout)
      use inland_parameters, only: lbeg, lend, numlv
      use inland_comwork, only: latscale, lonscale, latindex, lonindex
      use inland_lsmmapib, only: lati, loni
      implicit none

! Subroutine parameters
      integer imonthp, loopi, kpti, kptj, ipointout
      logical rdlsf

! Local variables
      integer i ! iterator variable :)

!   lati and loni variables: here we set lati/loni, both vectors npoi-sized,
! with the corresponding latscale and lonscale values for the given point.
!   For this to be done, we have to find a way to guess the corresponding
! latscale and lonscale to the given npoi.
      do 100 i = lbeg,lend
         lati(i)=latscale(latindex(i))
         loni(i)=lonscale(lonindex(i))
100   continue

!    imonthp: This is a new resource. INLAND did not support speci-
! fying restart month-wise, but only year-wise. So, for our approach, restart
! month is forced to be january (month 1).
      imonthp = 12 ! It denotes the last 'successfully ran' month.

!    kpti and kptj. While we are not parallel, kpti and kptj will be just the
! entire vector, in other words, from 1 to npoi.
      kpti = lbeg
      kptj = lend

!    Please note: lbeg is the first kpti of all intervals. lend is the last kptj
! of all intervals. If you have two intervals from 1 to 1000, you will have:
! 1st interval: kpti=1, kptj=500, lbeg=1, lend=1000
! 2nd interval: kpti=501,kptj=1000, lbeg=1, lend=1000
!    So, kpti/kptj is to sweep the particular set; lbeg/lend is to allocate
! the entire grid!

!    loopi indicates the current 'run' for the same process. Each process runs
! 'numlv' times, and its 'npoi' vector, split. If you have 4 processes running
! 8 times each, you will have the 1..npoi vector split by 4 intervals; each
! interval is split into 8 intervals which can be run in parallel. Each inter-
! val of these 8 is recognized by the variables 'begpnt' and advances up to 
! 'numpnt'. If this is not the last interval, numpnt+1 is the 'begpnt' for the 
! next interval.
!    If running a serial program, you have a single interval which is ran 1 time
! so no split is made at all and we have just a 1..npoi vector.

!    rdlsf is read from 'atm.parm' among with several other vars.
! It determines whether to initialize data for Antarctica and Groeland with spe-
! cial values, used only on coldstart subroutine.
      rdlsf = .true.

!    ipointout is read from 'atm.parm' among with several other variables.
! while we just need to reproduce SAGE's behavior, we will not read it but 
! assume it's default value
      ipointout = 0
end subroutine map_offline
