#include "inland_config.h"
#include "inland_compar.h"
#include "inland_comage.h"
#include "inland_comtex.h"

module inland_parameters
      implicit none
      public
      save

      integer rddaydims  !Gabriel Abrahao: added in the namelist declaration in main_offline and used in rdday. 0 for reading daily data in 4d, 1 in 3d
      integer isparse  !Gabriel Abrahao: added in the namelist declaration in main_offline and used to choose IO subroutines 0: Normal 3d/4d inputs 1: All inputs are one-dimensional vectors represeting a sparse matrix
      integer irdcropparmaps  !Gabriel Abrahao: added in the namelist declaration in main_offline and set to 1 if yearly crop parameter maps are to be used (e.g. planting dates)

      integer myid       ! mpi process (same as iam in ccm)
      integer nproc      ! number of mpi processes (same as npes in ccm)
      integer mpicomw    ! global mpi communicator
      integer iermpi     ! mpi error flag

      integer nl         ! number of lines of input file (0D Agro)

! Beginning and ending indices and ranges for each mpi process
!  Variables ending in 1 refer to actual land mesh
!  Variables not ending in 1 count multiplicity of land points
      integer, dimension(:), allocatable :: npbeg, npend, npnum
      integer, dimension(:), allocatable :: npbeg1, npend1, npnum1

! Beginning and ending indices for latitudinal decomp.
!  Variables ending in 1 refer to actual land mesh
!  Variables not ending in 1 count multiplicity of land points
      integer, dimension(:), allocatable :: njbeg, njend
      integer, dimension(:), allocatable :: njbeg1, njend1

! Land mesh index with multiplicity

      integer lbeg       ! local beginning meshpoint index
      integer lend       ! local ending meshpoint index
      integer lnum       ! local number of meshpoints

! Singular land mesh index

      integer lbeg1      ! local beginning meshpoint index
      integer lend1      ! local ending meshpoint index
      integer lnum1      ! local number of meshpoints

! define land surface 2-d grid and 1-d vector length

      integer mlpt     ! multiplicity of land points
      integer npoi1    ! number of land points on lsmlon x lsmlat grid
      integer npoi     ! number of land points (counting multiplicity)
! npoi used to be hard-coded (LPT in include/inland_compar.h), 
! but is now calculated in readit

! Grid resolution: longitude and latitude respective resolution in degrees
      real, parameter :: xres=XRES, yres=YRES

! Parallelization parameters
      integer numlv   !number of calls to inland (per process)
      integer mpt     !maximum number of points per call to inland (per process)
!! numlv is now set in inland_main_offline
!      parameter (numlv=1) ! not used on offline model as of yet

! define miscellaneous inland parameters
      integer nband   !number of solar radiation bands: vis, nir
      integer nsoilay !number of soil layers
      integer nsnolay !number of snow layers 
      integer npft    !number of plant functional types
      integer scpft   !starting index for crop pfts - C. Kucharik 
      integer ecpft   !ending index for crop pfts 
      integer npftu   !number of plant functional types in upper canopy
      integer ndat    !number of textural classes
      integer nvegtype!number of vegetation types - read in params/vegetation
      integer nlandusetype!number of land use types - read in params/vegetation

      parameter (nband=NUMBANDS, nsnolay=NUMSNOWLAYERS, &
                 npft=NUMPFT, scpft=SCPFT, ecpft=ECPFT, npftu=NUMPFT_UPPERCANOPY, ndat=NUMSOILTYPES)

! copies of atm variables
      integer, pointer :: plona, plata

! equivalents on inland model from SAGE

! longitude and latitude dimension of domain
      integer, parameter :: nlon=PLON, nlat=PLAT

      integer, target :: nlonsub, nlatsub
      integer ndaypm(12) ! Number of days per month
      integer ndaypy     ! Number of days per year

! number of days per month
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

! Former inland_control module (moved back to inland_parmeters which upgrades compar)

! some constants
      real*8, parameter :: pi=PISTOS

      real*8 :: epsilon, dtime, stef, vonk, grav, & !FIXME:dtime as floating pt?
                tmelt, hfus, hvap, hsub, ch2o,    &
                cice, cair, cvap, rair, rvap,     &
                cappa, rhow
!      epsilon        ! small quantity to avoid zero-divides and other
!                     ! truncation or machine-limit troubles with small
!                     ! values. should be slightly greater than o(1)
!                     ! machine precision.
!      dtime          ! inland time step (seconds)
!      stef           ! stefan-boltzmann constant (W m-2 K-4)
!      vonk           ! von karman constant (dimensionless)
!      grav           ! gravitational acceleration (m s-2)
!      tmelt          ! freezing point of water (K)
!      hfus           ! latent heat of fusion of water (J kg-1)
!      hvap           ! latent heat of vaporization of water (J kg-1)
!      hsub           ! latent heat of sublimation of ice (J kg-1)
!      ch2o           ! specific heat of liquid water (J deg-1 kg-1)
!      cice           ! specific heat of ice (J deg-1 kg-1)
!      cair           ! specific heat of dry air at constant pressure (J deg-1 kg-1)
!      cvap           ! specific heat of water vapor at constant pressure  (J deg-1 kg-1)
!      rair           ! gas constant for dry air (J deg-1 kg-1)
!      rvap           ! gas constant for water vapor  (J deg-1 kg-1)
!      cappa          ! rair/cair
!      rhow           ! density of liquid water (all types) (kg m-3)

! vzero(mpt): a real array of zeros, of length mpt 
      real*8, dimension(:), allocatable :: vzero
      real*8, dimension(:), allocatable :: garea

! model I/O variables, always access using trim()
! values can be overriden at runtime with env. variables 
! INLAND_INDIR, INLAND_INFILE, INLAND_DATADIR, INLAND_OUTDIR
!      character*255 indir, infile, datadir, outdir

! parameters identified on Optis efforts
! canopy
      real*8 stemae,     & ! stem activation energy
             rootae,     & ! root activation energy
             rrootpy,    & ! root respiration coefficient in kg/year units
             rwoodpy,    & ! wood respiration coefficient in kg/year units
             tempvm_nec, & ! tempvm's numerator's exponential coefficient
             rgrowth,    & ! gowth respiration coefficient (fraction)
             dispuhf,    & ! factor to weigh upper canopy parcel of z1 in dispu
             ialoglhf,   & ! initial alogl height fraction
             ialoguhf,   & ! initial alogu height fraction
             avmuir        ! average diffuse optical depth
! soil
      real*8 stressfac,  & ! to calculate moisture stress factor
             kfactor       ! mult. factor of decay constants for carb pools k**

end module inland_parameters
