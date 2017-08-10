#include "inland_config.h"
module inland_control
      implicit none
      public
      save

!
! Model control variables
!
      integer nrun      ! number of years in this simulation
      integer iyearout  ! 0: no yearly output 1: yearly output
      integer imonthout ! 0: no monthly output 1: monthly output
      integer idailyout   ! 0: no daily output 1: daily output
      integer ipointout ! 0: no point output
                        ! 1: point output
      integer isimveg   ! 0: static veg 
!                       ! 1: dynam veg initialized w/ fixed
!                       ! 2: dynam veg initialized w/ cold start
      integer isimfire  ! 0: fixed fire (0.5% prob.), 1: dynamic fire, 2: CTEM fire method, 3: no fire (default 0, ignored if isimveg=0)
      integer isimland  ! 0: fixed land  1: dynamic land
      integer isimco2   ! 0: fixed co2   1: changing co2
      integer ccmexist  ! 0: calculates existence arrays using climatological temperatures
!                       ! 1: calculates existence arrays using modeled monthly temperatures 
!
!Castanho Kai included isinfilt, isimrwu
      integer isinfilt  ! Infiltration Function
!                       ! 0: according to Darcy (1856); 
!                       ! 1: according to Green-Ampt (1911)
      integer isimrwu   !Root water uptake module
!                       ! 0: according to Foley et al., 1996; 
!                       ! 1 according to Li et al. (2006)
!!
      integer soicspin  ! 0: no spinup procedure for soil carbon
!                       ! 1: acceleration procedure used
!
!
      integer spinmax   ! degree of soil carbon acceleration
!
      real*8 spincons   ! number of times soilbgc subroutine is
!                       ! called for each day in simulation during
!                       ! acceleration procedure
      integer eqyears   ! years used past spin up to slowly bring
!                       ! entire soil to a steady equilibrium
      real*8 spinfrac   ! fraction of nspinsoil time used to
!                       ! spin up soil at maximum spin up rate
!
      integer nspinsoil ! year of simulation that soil c reaches equilibrium
!
      integer istep,  & ! timestep counter (per day)
              iyear0, & ! 1st year of simulation
              imonth0,& ! 1st month of simulation
              iday0,  & ! 1st day of simulation
              iyear,  & ! current year of simulation
              imonth, & ! current month of simulation
              iday,   & ! current day of simulation
              jday,   &  ! julian day of the simulation
              overveg
!
! Castanho HP included flags in ibisinfile to read Heterogeneous Parameterization

      integer itauw,	& !0:reads 1 dimension tauwood0 from parms canopy; 1:reads 2 dimension from input maps
              ivmax,	& !0:reads 1 dimension vmax_pft from parms canopy; 1:reads 2 dimension from input maps
              isla,	& !0:reads 1 dimension specla from parms canopy; 1:reads 2 dimension from input maps
              ica     	  !0:reads 1 dimension carbon allocation to wood, leaf and root from parms canopy; 1:reads 2 dimension from input maps


! model I/O variables, always access using trim()
! values can be overriden at runtime with env. variables 
! INLAND_INDIR, INLAND_INFILE, INLAND_DATADIR, INLAND_OUTDIR
      character*255 indir, infile, datadir, outdir

! ran2val optional, force ran2() return value for debug, default=0 to disable
! value can be overriden at runtime with env. variable INLAND_RANDOMVAL
      real*8 :: env_ran2val

! debug option to run the model quicker, run 1 day and 1 month only
! setup with env. variablt INLAND_FASTEXEC
      integer env_fastexec

! write real variables with single precision(float) - default is double precision
      integer env_floatout

! write compressed netcdf files (netcdf-4 with DEFLATE compression) - requires hdf-5
!   0: use uncompressed netcdf-classic format 
!   >1: use netcdf4-format with DEFLATE compression (requires hdf-5) 
!   value specifies compression level (1-9) (default 0, 2 is recommended)
!   note that compression is not recommended for small domains
      integer env_compressout

! define compressed file chunking
!   0: use netcdf default
!   1: lines (widthx1)
      integer env_chunkout

! debug output - set to higher than 0 for debug output, higher values increase verbosity
      integer env_debug 

! name of the vegtype file to read (and tileprop if mlpt>1)
! or default if blank (read in readit subroutine)
      character*1024 vegtypefile

! name of the hr map file to read, or none if blank (read in rdhrmapfile subroutine)
      character*1024 hrmapfile

! name of the crops map file to read, or none if blank (read in readit subroutine)
! TODO - should this be in comcrop?
      character*1024 cropsfile
! Variables required by land-use (rdlndtrans) subroutine

      integer iyrluc    ! 1st year to read land transition data
      integer nluc      ! number of years to read land transition data

#ifdef SINGLE_POINT_MODEL
      real nan          ! not a number representation (used on 0D)
      parameter(nan=-999.99)
#endif /* SINGLE_POINT_MODEL */
!
end module inland_control
