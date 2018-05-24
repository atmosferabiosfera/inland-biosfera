#include "inland_config.h"
!----------------------------------------------------------
logical function dir_exists(dir_name)

      implicit none

      character*255 dir_name
      logical :: dir_e

#if defined(__GFORTRAN__)
! gfortran does not support directory 
      inquire( file=trim(dir_name)//'/', exist=dir_e )
#elif defined(PGF90)
      inquire( file=trim(dir_name)//'/', exist=dir_e )
#else
! this works with ifort, test with others
      inquire( directory=trim(dir_name), exist=dir_e )
#endif
      dir_exists = dir_e

end function dir_exists


!----------------------------------------------------------
subroutine rd_env
!----------------------------------------------------------
!  Read various environment variables
! see docs/README.advanced for env. var. definition
!
!---------------------------Code history--------------------------------
! Author: Etienne Tourigny
! 
!----------------------------------------------------------

      use inland_parameters
      use inland_control

      implicit none
  
! model I/O variables
      character*255 chtmp
      logical :: dir_e
      logical dir_exists

#ifndef SINGLE_POINT_MODEL
#define INLAND_SUBDIR "offline/grid"
#define INLAND_INFILE "conf/inland-grid.infile"
#else
#define INLAND_SUBDIR "offline/single_point"
#define INLAND_INFILE "conf/inland-single_point.infile"
#endif

! Default path definitions, can be overriden by user at runtime
! FIXME is this ok in GNU only?
! INLAND_INDIR=data
! INLAND_INFILE=$INLAND_INDIR/<INLAND_INFILE>
! INLAND_DATADIR=./input or $INLAND_INDIR/$INLAND_SUBDIR/data
! INLAND_OUTDIR=.

      indir = './data/'//INLAND_SUBDIR
      datadir = './input'
      call getenv("INLAND_INDIR", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         indir = trim(chtmp)//'/'//INLAND_SUBDIR
         datadir = trim(indir)//'/data'
      end if
!      inquire( file=trim(indir)//'/', exist=dir_e )
      dir_e = dir_exists(indir)
      if ( .not. dir_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input directory '//trim(indir)//' does not exist!'
         write (*,*) 'make sure INLAND_INDIR is set to proper path'
         stop 1
      end if

      infile = trim(indir)//'/'//INLAND_INFILE
      call getenv("INLAND_INFILE", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         infile = trim(chtmp)
      end if
      inquire( file=trim(infile), exist=dir_e )
      if ( .not. dir_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input file '//trim(infile)//' does not exist!'
         write (*,*) 'make sure INLAND_INDIR or INLAND_INFILE is set to proper path (or not set)'
         stop 1
      end if

      call getenv("INLAND_DATADIR", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         datadir = chtmp
      end if
!      inquire( file=trim(datadir)//"/", exist=dir_e )
      dir_e = dir_exists(datadir)
!      inquire( directory=trim(datadir), exist=dir_e )
      if ( .not. dir_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: data directory '//trim(datadir)//' does not exist!'
         write (*,*) 'make sure INLAND_DATADIR is set to proper path'
         stop 1
      end if

      outdir = 'output'
      call getenv("INLAND_OUTDIR", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         outdir = chtmp
      end if
!      inquire( file=trim(outdir)//"/", exist=dir_e )
      dir_e = dir_exists(outdir)
      if ( .not. dir_e ) then
         call system( 'mkdir -p '//trim(outdir) )
      end if
!      inquire( file=trim(outdir)//"/", exist=dir_e )
      dir_e = dir_exists(outdir)
      if ( .not. dir_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: output directory '//trim(outdir)//' does not exist!'
         write (*,*) 'make sure INLAND_OUTDIR is set to proper path and directory exists'
         stop 1
      end if

#define CLAMPVAL(val,minval,maxval) \
      if ( val .lt. minval ) val = minval ; if ( val .gt. maxval ) val = maxval

! env_ran2val - set value from INLAND_RANDOMVAL if exists
      env_ran2val = 0
      call getenv("INLAND_RANDOMVAL", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         read (chtmp,*) env_ran2val
      end if
      CLAMPVAL( env_ran2val, 0, 1 )

! env_fastexec - set value from INLAND_FASTEXEC if exists
      env_fastexec = 0
      call getenv("INLAND_FASTEXEC", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         read (chtmp,*) env_fastexec
      end if
      if ( env_fastexec .lt. 0 ) env_fastexec = 0

! env_floatout - write real variables with single(float) precision (single by default)
      env_floatout = 1
      call getenv("INLAND_SINGLE_PRECISION", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         read (chtmp,*) env_floatout
      end if
      CLAMPVAL( env_floatout, 0, 1 )

! env_compressout - write compressed netcdf files (0-9)
      env_compressout = 0
      call getenv("INLAND_COMPRESSOUT", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         read (chtmp,*) env_compressout
      end if
      CLAMPVAL( env_compressout, 0, 9 )

! env_chunkout - define compressed variable chunking
      env_chunkout = 0
      call getenv("INLAND_CHUNKOUT", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         read (chtmp,*) env_chunkout
      end if
      CLAMPVAL( env_chunkout, 0, 1 )

! env_debug - set to higher than 0 for debug output
      env_debug = 0
      call getenv("INLAND_DEBUG", chtmp)
      if ( trim(chtmp) .ne. "" ) then
         read (chtmp,*) env_debug
      end if

      if ( env_debug .gt. 2 ) then
         write (*,9000) env_ran2val,env_fastexec,env_floatout,env_compressout,env_chunkout,env_debug
      end if
      9000  format (1x,'INFO: env vars                            : ',f8.2,5i2)
     
end subroutine rd_env
!*******************************************************************************
