#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine readit(isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
! ---------------------------------------------------------------------
! reads in initialization files and initializes some fields
      use inland_parameters
      use inland_combcs
      use inland_comsoi
      use inland_comveg
      use inland_compft                  ! Castanho HP, 2013
      use inland_comwork
      use inland_subgrid
      use inland_control, only: cropsfile, datadir, vegtypefile, hrmapfile, env_debug, itauw, ivmax, isla, ica  !Castanho HP, 2013
      use inland_comcrop
      use inland_comnitr

      implicit none

! Arguments
      integer isimveg, & ! dynamic vegetation (1) or static (0)
              iwest, jnorth

      real*8 snorth, ssouth, swest, seast

! Local variables
      integer  istat,  & ! netcdf error flag
           i, j, k, ntime,& ! loop indices
           jj, ii,     &
           ndim,       & ! number of dimensions
           nlpoints      ! number of land points

      integer tmpint
      integer, dimension (1) :: max_loc

      real*8 xlat, totreal
      real*8 caverf ! Castanho HP, 2013 for verification of awood+aroot+aleaf=1

      integer jsouth, ieast
      integer ilpt, ipt, ipt2  ! tmp point indices for mlpt loop
      integer mlpt1 !mlpt used to read vegtype
      real*8, dimension(:,:), allocatable :: buffer2 ! used to read cropsfile

! xmask(nlon,nlat) ::> equivalence --> pointer 
! cant use equivalence ( xmask(1,1), work(1) ) to xmask with allocatable arrays
!   The work variable is ndim2 which is nlon*nlat size. Here, we just declare
! xmask as (nlon,nlat)-dimensioned and use it on a new memory space.
!   If there's still problem with that approach, we will have to allocate work 
! statically. As of yet, nlon and nlat are compile-time options and this
! will not be a problem. Even if we were to get nlon and nlat from input config
! file (inland-grid.infile/inland-single_point.infile), this would already have been 
! done.
      real*8 xmask(nlon,nlat)

      character*1024 filen, filen2
      integer istart(4), icount(4) ! for reading vars
      logical :: file_e


! ---------------------------------------------------------------------
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /

! 2-d surface and vegetation arrays

#ifndef SINGLE_POINT_MODEL
! land mask, latitudes, and longitudes

      filen = trim(datadir)//'/surta.nc'     

      ! make sure this file exists, if not print error and exit
      inquire( file=trim(filen), exist=file_e )
      if ( .not. file_e ) then
         write (*,*) ''
         write (*,*) 'ERROR: input file '//trim(filen)//' does not exist!'
         write (*,*) 'make sure INLAND_DATADIR is set to proper path and that file exists'
         stop 1
      end if

      ! first read lon/lat vars directly
      aname = 'longitude'
      icount = (/nlon,1,1,1/)
      call readvar(filen,aname,lonscale,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999

      aname = 'latitude'
      icount = (/nlat,1,1,1/)
      call readvar(filen,aname,latscale,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999

      ! now read other vars, reset icount
      icount = (/nlon,nlat,1,1/)

      aname = 'surta'
      call readvar(filen,aname,xmask,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999

      if (abs(abs(lonscale(1)-lonscale(2))-xres).gt.0.001 .or. abs(abs(latscale(1)-latscale(2))-yres).gt.0.001) then
         write (*,9000)
         write(*,*) 'resolution mismatch!'
         write(*,*) 'xres, yres in inland_compar.h = ', xres, yres
         write(*,*) 'xres, yres from '//trim(datadir)//'/surta.nc = ', &
                    abs(lonscale(1)-lonscale(2)), abs(latscale(1)-latscale(2))
         stop 1
      end if
#endif /* SINGLE_POINT_MODEL */

! subset the grid if not default whole grid
      if (snorth.lt.latscale(1) .or. ssouth.gt.latscale(nlat) .or. &
          swest.gt.lonscale(1) .or.  seast.lt.lonscale(nlon)) then
         jnorth = 0
        
         if (snorth .lt. (latscale(nlat)+latscale(nlat-1))/2.) then
            jnorth = nlat
         elseif (snorth .ge. (latscale(1)+latscale(2))/2.) then
            jnorth = 1
         else
            do j = nlat-1,1,-1
               if (snorth .ge. (latscale(j)+latscale(j+1))/2.) jnorth = j
            end do
         end if
         jsouth = 0
         if (ssouth .lt. (latscale(nlat)+latscale(nlat-1))/2.) then
            jsouth = nlat
         elseif (ssouth .ge. (latscale(1)+latscale(2))/2.) then
            jsouth = 1
         else
            do j = nlat-1,1,-1
               if (ssouth .ge. (latscale(j)+latscale(j+1))/2.) jsouth = j
            end do
         end if
!
         iwest = 0
         if (swest .lt. (lonscale(1)+lonscale(2))/2.) then
            iwest = 1
         elseif (swest .ge. (lonscale(nlon)+lonscale(nlon-1))/2.) then
            iwest = nlon
         else
            do i = 2, nlon
               if(swest .ge. (lonscale(i)+lonscale(i-1))/2.) iwest=i
            end do
         end if
         ieast = 0
         if (seast .lt. (lonscale(1)+lonscale(2))/2.) then
            ieast = 1
         elseif (seast .ge. (lonscale(nlon)+lonscale(nlon-1))/2.) then
            ieast = nlon
         else
            do i = 2, nlon
               if(seast .ge. (lonscale(i)+lonscale(i-1))/2.) ieast=i
            end do
         end if
         nlonsub = ieast - iwest + 1
         nlatsub = jsouth - jnorth + 1
         istart(1) = iwest
         icount(1) = nlonsub
         istart(2) = jnorth
         icount(2) = nlatsub
      else
         iwest = 1
         ieast = nlon
         jnorth = 1
         jsouth = nlat
         nlonsub = nlon
         nlatsub = nlat
      end if

! Variables in inland parameter used for readit subroutine.
! Attach (via pointers) nlonsub/nlatsub to plona/plata from inland_compar
      plona => nlonsub
      plata => nlatsub
! The lmask(nlonsub,nlatsub) or rather lmask(plon,plat) must be allocated here!
! Remember plona => nlonsub and plata => nlatsub on inland_prealloc!
      allocate(lmask(plona,plata))
      lmask(:,:) = 0

! First, we calculate nlpoints
! At this point we cycle thru all the loop without changing anything.
! As in sage, here i/j refers to entire grid (1 to nlon/nlat), 
! ii/jj refers to subgrid (1 to nlonsub/nlatsub)
      nlpoints = 0
      do j = jnorth, jsouth
         jj = j - jnorth + 1
         do i = iwest, ieast
            ii = i - iwest + 1
#ifdef SINGLE_POINT_MODEL
            if ((ii*jj).eq.1) then
               lmask(ii,jj) = 1
            else
               lmask(ii,jj) = 0
            endif
#else /* SINGLE_POINT_MODEL */
            lmask(ii,jj) = nint(xmask(i,j))
#endif /* SINGLE_POINT_MODEL */
            if (lmask(ii,jj).eq.1) then
               nlpoints = nlpoints + 1
            end if
         end do !i
      end do !j

      write (*,9010) 
      write (*,9020) nlpoints
      if ( mlpt .gt. 1 ) write (*,9030) mlpt
      write (*,9010) 

! now we set npoi, lbeg/lend and mpt here, this should change with openmp/mpi
! this used to be hard-coded in include/inland_compar.h (LPT variable)
      npoi1 = nlpoints
      npoi = mlpt * npoi1
      lbeg = 1
      lend = npoi
      mpt = lend - lbeg + 1

      if ( lend .ne. npoi ) then
         write(*,*) 'ERROR in readit - lend != npoi - must rewrite code that deals with npoi'
         stop 1
      end if

! ---------------------------------------------------------------------
! allocate variables dependant on npoi which will be used here, was in inland_prealloc
! TODO - should this go into a new subroutine ?
! or should we add to inland_alloc and call from here ?

! Allocate variables in comwork
!--------------------------------
      allocate(lonindex(lbeg:lend),latindex(lbeg:lend))
      lonindex(:) = 0
      latindex(:) = 0
      allocate(buffer(lbeg:lend))
      allocate(buffer1(1:npoi1))

! Allocate variables in comveg
!-------------------------------
      allocate(landusetype(lbeg:lend))
      allocate(croptype(lbeg:lend))
      landusetype(:) = 1.
      croptype(:) = 0.

! Allocate variables in subgrid
!-------------------------------
! subgrid tiles - this has to change if lbeg != 1 and lend != npoi
! defaults are tilefrac=1 for dominant tile, 0 for others
      allocate(tilefrac(1:npoi))
      tilefrac(:) = 0.
      allocate(waterfrac(1:npoi1))
      waterfrac(:) = 0.
      maxsubgrid = mlpt !should be user-config, but set to mlpt for now
      allocate(itilechild(1:npoi1,maxsubgrid))
      itilechild(:,:) = 0
      allocate(ntilechild(1:npoi1))
      ntilechild(:) = 0
      allocate(itileparent(1:npoi))
      itileparent(:) = 0
      ! first subgrid is 1:npoi1, this CANNOT change because some code depends on this
      do i = 1,npoi1
         tilefrac(i) = 1.
         ipt2 = subgrid_set_index(i,1,i)
         if (ipt2 .eq. 0) then
            write(*,*) 'ERROR subgrid_set_index(',i,1,i,') returned ',ipt2
         end if
      end do

! Allocate variables in combcs
!-------------------------------
!   Please note: lmask have to be allocated into readit as its dimention is
! defined just before lmask is used!
      allocate(xintopo(lbeg:lend),xinveg(lbeg:lend),deltat(lbeg:lend), &
               xint(lbeg:lend,12))
      xintopo(:) = 0.
      xinveg(:) = 0.
      deltat(:) = 0.
      xint(:,:) = 0.

! combcs's readit-specific variables
      allocate(clmwet(lbeg:lend,12),clmtrng(lbeg:lend,12))
      allocate(obswet(lbeg:lend,12),obstrng(lbeg:lend,12))

! combcs's anomaly and readit specific variables
      allocate(clmt(lbeg:lend,12),clmprec(lbeg:lend,12),clmcld(lbeg:lend,12), &
               clmq(lbeg:lend,12),clmwind(lbeg:lend,12))
      clmt(:,:)=0.
      clmprec(:,:)=0.
      clmcld(:,:)=0.
      clmq(:,:)=0.
      clmwind(:,:)=0.

      allocate(obst(lbeg:lend,12),obsprec(lbeg:lend,12),obscld(lbeg:lend,12), &
               obsq(lbeg:lend,12),obswind(lbeg:lend,12))
      obst(:,:)=0.
      obsprec(:,:)=0.
      obscld(:,:)=0.
      obsq(:,:)=0.
      obswind(:,:)=0.

! combcs's weather generator and readit variables
! Please note these are read at readit thus no need to assign zeroes.
      allocate(xinwind(lbeg:lend,12),xintrng(lbeg:lend,12), &
               xinprec(lbeg:lend,12),xincld(lbeg:lend,12),xinq(lbeg:lend,12), &
               xinwet(lbeg:lend,12))

      allocate(xinwindmon(lbeg:lend,12),xintrngmon(lbeg:lend,12), &
               xinprecmon(lbeg:lend,12),xincldmon(lbeg:lend,12),xinqmon(lbeg:lend,12), &
               xinwetmon(lbeg:lend,12),xintmon(lbeg:lend,12))

! Castanho HP, 2013 combcs's included new input heterogeneous parameterization maps read from input folder
!------------------------------------------------------------------------------------------------------
      allocate(xintauwood(lbeg:lend,npft), xinawood(lbeg:lend,npft), xinaroot(lbeg:lend,npft), &
                 xinvmax(lbeg:lend,npft), xinspecla(lbeg:lend,npft),xinaleaf(lbeg:lend,npft), &
                 vmax_pft(lbeg:lend,npft), tauwood0(lbeg:lend,npft), awood(lbeg:lend,npft),    &
                 aroot(lbeg:lend,npft),astem(lbeg:lend,npft),aleaf(lbeg:lend,npft),specla(lbeg:lend,npft))       
      tauwood0(:,:) = 0.
      vmax_pft(:,:) = 0.     
      awood(:,:) = 0. 
      aroot(:,:) = 0. 
      astem(:,:) = 0. 
      aleaf(:,:) = 0. 
      specla(:,:) = 0.     

! Variables in comsoi used for readit subroutine
      allocate(sand(lbeg:lend,nsoilay),clay(lbeg:lend,nsoilay))
      sand(:,:) = 0.
      clay(:,:) = 0.

      allocate(garea(lbeg:lend))
      garea(:) = 0.

! gabriel abrahao: Allocate the crop parameters that can also be read as maps (rdcropparmaps) here, and replicate their temporary values to all points. If rdcropparmaps is called later, those will be overriden
      if (isimagro.gt.0) then
         allocate(pmmin(lbeg:lend,npft),pdmin(lbeg:lend,npft))
         do j = scpft, ecpft 
            pdmin(:,j) = pdmin_temp(j)
            pmmin(:,j) = pmmin_temp(j)
         end do
      end if

! ---------------------------------------------------------------------


! Now we write data to the structures.
! initialize lonindex, latindex for use in arr2vec, vec2arr, etc.
! and calculate the approximate the area of the land gridcells
! here, i/j refers to entire grid (1 to nlon/nlat), 
! ii/jj refers to subgrid (1 to nlonsub/nlatsub)
      nlpoints = 0
      do j = jnorth, jsouth
         jj = j - jnorth + 1
         do i = iwest, ieast
            ii = i - iwest + 1
#ifdef SINGLE_POINT_MODEL
            if ((ii*jj).eq.1) then
               lmask(ii,jj) = 1
            else
               lmask(ii,jj) = 0
            endif
#else /* SINGLE_POINT_MODEL */
            lmask(ii,jj) = nint(xmask(i,j))
#endif /* SINGLE_POINT_MODEL */
            if (lmask(ii,jj).eq.1) then
               nlpoints = nlpoints + 1
               lonindex(nlpoints) = ii
               latindex(nlpoints) = jj
               xlat = latscale(j) * pi / 180.0
               garea(nlpoints) = yres * 111400.0 * xres * 111400.0 * cos(xlat)

               ! create subgrid indexes
               ! replicate lonindex, latindex and garea to subgrid tiles
               do ilpt= 1,mlpt
                  ipt = nlpoints+(ilpt-1)*npoi/mlpt
                  ! set subgrid indexes
                  ipt2 = subgrid_set_index(nlpoints,ilpt,ipt)
                  if (ipt2 .eq. 0) then
                     write(*,*) 'ERROR subgrid_set_index(',nlpoints,ilpt,ipt,') returned ',ipt2
                  else if (ipt .ne. ipt2) then
                     write(*,*) 'ERROR subgrid_set_index(',nlpoints,ilpt,ipt,') returned',ipt2, 'should be ',ipt
                  end if
                  lonindex(ipt) = ii
                  latindex(ipt) = jj
                  !FIXME calculate based on actual fraction, for now copy parent's area
                  garea(ipt) = garea(nlpoints)
               end do

            end if
         end do !i
      end do !j

#ifdef SINGLE_POINT_MODEL
! NOTE NOTE NOTE ::: Single point model will return from here. Nothing below
!  this point will be executed!

      do j = 1, npft
         do i = 1,npoi
            tauwood0(i,j)= tauwood0p(j)
            specla(i,j)= speclap(j)
            awood(i,j)= awoodp(j)
            aroot(i,j)= arootp(j)
            aleaf(i,j)= aleafp(j)
            vmax_pft(i,j)= vmax_pftp(j)
         enddo
      enddo

      return

#else /* SINGLE_POINT_MODEL */
      do j = jnorth, jsouth
         jj = j - jnorth + 1
         latscale(jj) = latscale(j)
      end do
      do i = iwest, ieast
         ii = i - iwest + 1
         lonscale(ii) = lonscale(i)
      end do

! ---------------------------------------------------------------------
! now we read topo, vegtype and clim input

! dummy variable example, 4-d, but 3rd dim ('level') = 1
! copy and chanve for a new variable
!
!      filen = 'input/dummyv.nc'
!      aname = 'dummyv'
!      call readvar(filen,aname,'level',istart,icount,
!     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
!      if (istat.lt.0) then
!         write(*,9000)
!         print *, 'while reading dummyv'
!         stop 1
!      end if
!      call arr2vec(cdummy, xindummy)
!
! topography
      filen = trim(datadir)//'/topo.nc'
      aname = 'topo'

! note that this comment no longer applies, as readvar does not set dimension values
!   Remember ndim4 is the biggest between nlon,nlat,nband,nsoilay,nsnolay,npft;
! ndim2 (size of work) is nlon*nlat so as long as nlon and nlat are greater than
! or equal to 3, there will be no problem on the usage below. ndim2 is nlat*nlon

      call readvar(filen,aname,xintopo,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! fixed vegetation map
      if (isimveg .le. 1) then

         ! read in vegtypes, using tiles if mlpt > 1
         ! for information on the creation of vegtype files with tiles see docs/README.advanced

         mlpt1 = mlpt
         if ( trim(vegtypefile) .ne. '' ) then
            filen = trim(datadir)//'/'//trim(vegtypefile)
         else
            ! if mlpt == 1 or we read crops map file (isimagro==2), don't read vegtype+tileprop from tile files
            if ( ( mlpt == 1 ) .or. ( isimagro .eq. 2 ) ) then
               mlpt1 = 1
               filen = trim(datadir)//'/vegtype.nc'
            else
               if ( mlpt <= 9 ) then
                  write(filen, '(A,i1,A)'), trim(datadir)//'/tiles/vegtype-mlpt0',mlpt,'.nc'
               else
                  write(filen, '(A,i2,A)'), trim(datadir)//'/tiles/vegtype-mlpt',mlpt,'.nc'
               end if
               ! make sure file exists
               inquire( file=trim(filen), exist=file_e )
               if ( .not. file_e ) then
                  ! if this file does not exist, use normal file (tilefrac will be 1 only for first tile)
                  ! this is disabled by default, comment below to allow (used in dev. only!)
                  write (*,*) 'ERROR: vegtype input file '//trim(filen)//' does not exist!'
                  stop 1
                  mlpt1 = 1
                  filen = trim(datadir)//'/vegtype.nc'
               end if
            end if
         end if

         inquire( file=trim(filen), exist=file_e )
         if ( .not. file_e ) then
            write (*,*) ''
            write (*,*) 'ERROR: vegtype input file '//trim(filen)//' does not exist!'
            stop 1
         end if

         if ( env_debug .gt. 0 ) print *,'reading vegtype from '//trim(filen)

         ! icount(3) is the mlpt layers used in vegtype and tilefrac
         icount(3) = mlpt1
         icount(4) = 1

         aname = 'vegtype'
         call readvar(filen,aname,xinveg,istart,icount,mlpt,istat)
         if (istat.lt.0) goto 9999

         ! read subgrid tile fractions if mlpt > 1
         ! for creating/modyfing tile files see "Subgrid Tiling" section in docs/README.advanced
         if ( mlpt1 .gt. 1 ) then
            aname = 'tilefrac'
            call readvar(filen,aname,tilefrac,istart,icount,mlpt,istat)
            if (istat.lt.0) goto 9999

            ! if tilefrac contains mlpt+1 levels, last level is waterfrac
            if ( .true. ) then
               istart(3) = mlpt1+1
               icount(3) = 1
               call readvar(filen,aname,waterfrac,istart,icount,1,istat)
               ! TODO detect #levels, for now catch any error in readvar
               !if (istat.lt.0) goto 9999
               if (istat.lt.0) then
                  print *, 'warning: tilefrac in file '//trim(filen)//' should contain mlpt+1 level (last level is waterfrac)'
                  waterfrac(:) = 0.
               end if
               istart(3) = 1
               icount(3) = mlpt1
            end if
         end if ! tilefrac

         ! read HR map if requested
         ! TODO - don't do previous steps in this case - but for now use this as default
         if ( trim(hrmapfile) .ne. '' ) then
#ifdef HRMAP
            if ( env_debug .gt. 0 ) print *,'reading hrmap from '//trim(filen)
            call rdhrvegmap(trim(datadir)//'/'//trim(hrmapfile),'vegtype')
#else
            write(*,*) 'ERROR! hrmapfile defined in infile, but model not configured with --enable-highres_map option'
#endif
         end if

         ! set landuse type in function of xinveg
         do j = 1, npoi
            landusetype(j) = landusetypemap(nint(xinveg(itileparent(j))))
         end do

         ! agro
         ! when a point is defined as a specific crop, then set
         !   landusetype=2 (cropland) and croptype=<crop_pft> where crop_pft is the corresponding
         !   pft number, as defined in params/vegetation
         ! there are 2 ways to define crops, depending on isimagro value
         ! TODO add isimagro check in inland_test - an error happens if isimagro=-1
         if ( isimagro .gt. 0 ) then         

            ! if we using a unique crop (using icroptype)
            if ( isimagro .eq. 1 ) then

               if ( (icroptype .lt. scpft) .or. (icroptype .gt. ecpft) ) then
                  write(*,*) 'ERROR! icroptype bust be between',scpft,'and',ecpft
                  goto 9999
               end if

               ! make sure iwheattype is properly defined if using wheat
               !TODO put before
               if ( icroptype .eq. 15 ) then
                  if ( iwheattype .eq. 0 ) then
                     write(*,*) 'ERROR! iwheattype must be equal to 1 or 2'
                     goto 9999
                  end if
               end if

               if ( icroptype .eq. 13 .or. icroptype .eq. 14 .or. icroptype .eq. 16 ) then
                  if ( iwheattype .eq. 1 .or. iwheattype .eq. 2 ) then
                     write(*,*) 'ERROR! iwheattype must be equal to 0'
                     goto 9999
                  end if
               end if

               write (*,*) 'INFO: running in agro mode, all points have croptype:',icroptype

               xinveg(:) = 17 ! cropland vegtype
               landusetype(:) = 2
               croptype(:) = icroptype

            ! or read crops map file (cropsfile)
            else if ( isimagro .eq. 2 ) then

               if ( mlpt .eq. 1 ) then
                  write (*,*) 'ERROR! mlpt must be > 1 when reading cropsfile'
                  goto 9999
               end if

               if ( trim(cropsfile) .eq. '' ) then
                  write(*,*) 'ERROR! cropsfile must be defined if isimagro=2'
                  goto 9999
               end if

               filen = trim(datadir)//'/'//trim(cropsfile)
               aname = 'cropdata'
               
               write (*,*) 'INFO: running in agro mode, reading crops from '//trim(filen)
               
               ! number of agro pfts - cropsfile must have (at least) this number of crops defined 
               ! and in the same order as crop pfts 
               tmpint = ecpft - scpft + 1 
               allocate(buffer2(1:npoi1,tmpint))
               buffer2(:,:) = 0

               icount(3) = tmpint
               icount(4) = 1
               call readvar(filen,aname,cdummy,istart,icount,-1,istat)
               if (istat.lt.0) goto 9999
               do j = 1, tmpint
                  call arr2vec_tile(cdummy((j-1)*nlonsub*nlatsub + 1), buffer2(1,j), 1)
               end do

               ! specific case - one subgrid for each crop
               if ( mlpt .eq. (tmpint+1) ) then
                  do i = 1, npoi1
                     do j = 1, tmpint
                        ipt = subgrid_get_index(i,j+1)
                        ! first subgrid is natural vegetation, next subgrids are crops
                        xinveg(ipt) = 17 ! cropland vegtype
                        tilefrac(ipt) = buffer2(i,j)
                        tilefrac(i) = tilefrac(i) - tilefrac(ipt)
                        landusetype(ipt) = 2
                        croptype(ipt) = j + scpft - 1
                     end do ! tmpint
                  end do ! npoi1
                  
               ! general case - less subgrids than # of crops
               else
                  do i = 1, npoi1
                     ! first get all values so we can sort crops by fraction
                     do j = 1, tmpint
                        buffer(j) = buffer2(i,j) 
                     end do
                     ! now find maxes and use them as subgrids, in decreasing fraction
                     do j = 1, mlpt-1
                        ipt = subgrid_get_index(i,j+1)
                        max_loc = maxloc( buffer(1:tmpint) )
                        k = max_loc(1)
                        if ( buffer(k) .gt. 0. ) then
                           ! first subgrid is natural vegetation, next subgrids are crops
                           xinveg(ipt) = 17 ! cropland vegtype
                           tilefrac(ipt) = buffer(k)
                           tilefrac(i) = tilefrac(i) - tilefrac(ipt)
                           landusetype(ipt) = 2
                           croptype(ipt) = k + scpft - 1
                           buffer(k) = 0. ! remove from buffer
                        else
                           ! make this tile a natural, but tilefrac=0
                           xinveg(ipt) = xinveg(i)
                        end if
                     end do ! mlpt-1 
                  end do ! npoi1
               end if
               
               icount(3) = mlpt
               istart(3) = 1
               icount(4) = 1

               deallocate(buffer2)

            else

               write(*,*) 'ERROR! isimagro must be between 0 and 2'
               goto 9999

            end if ! cropsfile

         end if ! isimagro

         ! sanity checks on iwheattype
         if ( iwheattype .gt. 0 .and. icroptype .ne. 15 ) then
            !write(*,*) 'NOTICE: iwheattype set to 0 because icroptype != 15'
            iwheattype = 1
         end if
         if ( iwheattype .eq. 0 .and. icroptype .eq. 15 ) then
            write(*,*) 'NOTICE: iwheattype set to 1 because icroptype == 15'
            iwheattype = 1
         end if

         if ( mlpt .gt. 1 ) then

            ! make sure totals = 1.0 for each point
            ! verify with 
            ! cdo vertsum inland-tiles-1981.nc  tmp1.nc ;  cdo selvar,tilefrac tmp1.nc tmp2.nc ; cdo selvar,waterfrac inland-tiles-1981.nc tmp3.nc ; cdo add tmp2.nc tmp3.nc tmp4.nc
            do i = 1,npoi1
               totreal = dble(0)
               ! calculate grid total
               if ( env_debug .gt. 3 ) print *,''
               do ilpt = 1, mlpt
                  ipt = subgrid_get_index(i,ilpt)
                  if ( ipt .ne. 0 ) then
                     totreal = totreal + tilefrac(ipt)
                  end if
                  if ( env_debug .gt. 3 ) print *,i,ilpt,ipt,tilefrac(ipt),waterfrac(i),totreal
               end do
               totreal = totreal + waterfrac(i)
               ! if total > 1.0, exit
               if ( ( totreal - dble(1) ) .gt. 0.0001 ) then
                  write(*,9000)
                  write(*,9040),i,totreal
                  stop 1
               ! if total < 1.0, add missing to first tile (assuming first is dominant tile)
               ! TODO find dominant tile (i.e. the one with highest value)
               else if ( ( dble(1) - totreal ) .gt. 0.0001 ) then
                  ipt = subgrid_get_index(i,1)
                  if ( env_debug .gt. 0 ) write(*,9050),1,i,tilefrac(ipt),( tilefrac(ipt) + dble(1) - totreal )
                  ! if using hrmap, exit because this shouldn't happen 
                  ! and ihrtileparent indexes will be wrong
                  if ( trim(hrmapfile) .ne. '' ) then
                     write(*,9040),i,totreal
                     stop 2
                  end if
                  tilefrac(ipt) = tilefrac(ipt) + dble(1) - totreal
                  totreal = dble(1)
               ! adust for any rounding errors
               else if ( abs(totreal - dble(1) ) > 0.00000000000001 ) then
                  if ( env_debug .gt. 1 ) then
                     write(*,9060),i,totreal
                     write(*,9051),i,ipt,tilefrac(i),tilefrac(i) + dble(1) - totreal
                  end if
                  tilefrac(i) = dble(tilefrac(i)) + dble(1) - totreal
                  totreal = dble(1)
               end if
               if ( env_debug .gt. 3 ) print *,i,tilefrac(i),totreal,waterfrac(i)
            end do

            ! update ntilechild
            ntilechild(:) = 0
            do i = 1,npoi1
               tmpint = 0
               ! calculate total
               do ilpt = 1, mlpt
                  ipt = subgrid_get_index(i,ilpt)
                  if ( tilefrac(ipt) .gt. 0.0 ) then
                     tmpint = tmpint + 1
                  end if
               end do
               ntilechild(i) = tmpint
            end do

         end if ! mlpt > 1
         
      end if ! (isimveg .le. 1)          

      if ( env_debug .ge. 5 ) then
         print *,'xinveg:     ', xinveg
         if ( mlpt .gt. 1 ) print *,'tilefrac:   ', tilefrac
         print *,'landusetype:', landusetype
         print *,'croptype:   ' ,croptype
      end if

      ! make sure that xinveg != 0 for all points
!      do i = 1, npoi
!         if ( xinveg(i) .eq. 0 ) then
!            write (*,*) ''
!            write (*,*) 'ERROR: xinveg=0 for point',i
!            stop 1
!         end if
!      end do

!c 2-d soil array
!c
!     filen = 'input/soil.nc'
!     aname = 'soil'
!     call arr2vec(cdummy, soita)
!
! delta t
      icount(3) = 1
      icount(4) = 1
      filen = trim(datadir)//'/deltat.nc'
      aname = 'deltat'
      call readvar(filen,aname,deltat,istart,icount,0,istat)
      if (istat.lt.0) goto 9999

! 3-d soil texture array
!
! icount(3) is the 6 layers used in soita.sand.nc soita.clay.nc

      !gabriel abrahao this used to be 6
      if ( nsoilay .ne. 12 ) then
         write (*,*) ''
         write (*,*) 'WARNING: nsoilay = ',nsoilay,' but should be 12!'
         write (*,*) 'fix me in inland_readit'
         !stop 1
      end if

      !icount(3) = 6 !nsoilay 
      !gabriel abrahao this used to be 6
      icount(3) = 12 !nsoilay in input file
      icount(4) = 1
      filen = trim(datadir)//'/soita.sand.nc'
      aname = 'sandpct'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do j = 1, nsoilay
         call arr2vec(cdummy((j-1)*nlonsub*nlatsub + 1), sand(1,j))
      end do

      filen = trim(datadir)//'/soita.clay.nc'
      aname = 'claypct'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do j = 1, nsoilay
         call arr2vec(cdummy((j-1)*nlonsub*nlatsub + 1), clay(1,j))
      end do

! 3-d climate arrays
      icount(3) = 1
      icount(4) = 12

      filen = trim(datadir)//'/wetd.mon.nc'
      aname = 'wetd'

      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),clmwet(1,ntime))
      end do

      filen = trim(datadir)//'/temp.mon.nc'
      aname = 'temp'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),clmt(1,ntime))
      end do

      filen = trim(datadir)//'/trange.mon.nc'
      aname = 'trange'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),clmtrng(1,ntime))
      end do

      filen = trim(datadir)//'/prec.mon.nc'
      aname = 'prec'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),clmprec(1,ntime))
       end do

      filen = trim(datadir)//'/wspd.mon.nc'
      aname = 'wspd'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),xinwind(1,ntime))
      end do

      filen = trim(datadir)//'/cld.mon.nc'
      aname = 'cld'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),clmcld(1,ntime))
      end do

      filen = trim(datadir)//'/rh.mon.nc'
      aname = 'rh'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,12
         call arr2vec(cdummy((ntime-1)*nlonsub*nlatsub + 1),clmq(1,ntime))
      end do

! Castanho HP, 2013 reading from input the spatial heterogeneous map of aroot, awood, (aleaf = 1-awood-aroot) tauwood, specla and Vmax

! 3-d plant functional type array
! icount(3) is the 12 pft  used in residence time
      icount(3) = npft
      icount(4) = 1
      if (itauw.eq.0) then
      do j = 1, npft
         do i = 1,npoi
            tauwood0(i,j)= tauwood0p(j)
	 enddo
      enddo
      else
      filen = trim(datadir)//'/tauwood.nc'
      aname = 'tauwood'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do ntime = 1,npft
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1), xintauwood(1,ntime))
      enddo
      ndim = npoi*npft
      call scopya(ndim, xintauwood,tauwood0)
      end if

      if (ica.eq.0) then
      do j = 1, npft
         do i = 1,npoi
            awood(i,j)= awoodp(j)
            aroot(i,j)= arootp(j)
            aleaf(i,j)= aleafp(j)
         enddo
      enddo
      else
      filen = trim(datadir)//'/awood.nc'
      aname = 'awood'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do  ntime = 1,npft
        call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1), xinawood(1,ntime))
      enddo
      ndim = npoi*npft
      call scopya (ndim,xinawood,awood)

      filen = trim(datadir)//'/aroot.nc'
      aname = 'aroot'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do  ntime = 1,npft
        call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1), xinaroot(1,ntime))
      enddo
      ndim = npoi*npft
      call scopya (ndim,xinaroot,aroot)

      filen = trim(datadir)//'/aleaf.nc'
      aname = 'aleaf'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do  ntime = 1,npft
        call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1), xinaleaf(1,ntime))
      enddo
      ndim = npoi*npft
      call scopya (ndim,xinaleaf,aleaf)      
      end if
      
! Verifing if the sum of aleaf+awood+aroot=1 
!
   if(isimagro .eq. 0)then
      do i=1, npoi
         do j=1,12
            caverf=aleaf(i,j)+ awood(i,j)+aroot(i,j)
            if(caverf.gt.1.or.caverf.lt.1) then 
               print*,'Unexpected values : The total carbon allocation aroot + awood + aleaf is greater or lower then one' 
               print*,'ERROR in subroutine readit'
               stop
            endif 
         enddo
      enddo
    endif

      if (ivmax.eq.0) then
         do j = 1, npft
            do i = 1,npoi
               vmax_pft(i,j)= vmax_pftp(j)
            enddo
         enddo
      else
      filen = trim(datadir)//'/vmax.nc'
      aname = 'vmax'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do  ntime = 1,npft
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1), xinvmax(1,ntime))
      enddo
      ndim = npoi*npft
      call scopya (ndim,xinvmax,vmax_pft)
      endif

      if (isla.eq.0) then
         do j = 1, npft
            do i = 1,npoi
               specla(i,j)= speclap(j)
            enddo
         enddo
      else

      filen = trim(datadir)//'/specla.nc'
      aname = 'specla'
      call readvar(filen,aname,cdummy,istart,icount,-1,istat)
      if (istat.lt.0) goto 9999
      do  ntime = 1,npft
        call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1), xinspecla(1,ntime))
      enddo
      ndim = npoi*npft
      call scopya(ndim, xinspecla,specla)
      end if

! copy all 5 climatology fields to clm+anom fields for spin up
      ndim = npoi*12
      call scopya(ndim, clmt, xint)
      call scopya(ndim, clmtrng, xintrng)
      call scopya(ndim, clmprec, xinprec)
      call scopya(ndim, clmcld, xincld)
      call scopya(ndim, clmq, xinq)
      call scopya(ndim, clmwet, xinwet)
      
! return to main program
      return
#endif /* SINGLE_POINT_MODEL */

9000  format (1x,'ERROR in subroutine readit')
9010  format (1x,' ')
9020  format (1x,'INFO: number of land points               : ', i8)
9030  format (1x,'INFO: mult.  of land points (tiles)       : ', i8)
9040  format (1x,'ERROR: tilefrac total of point ',i8,' is ',f8.4,' (expecting 1.0)')
9050  format (1x,'INFO: tilefrac of tile ',i8,' of point ',i8,' changed from ',f6.4,' to ',f6.4)
9051  format (1x,'INFO: tilefrac of tile ',i8,' of point ',i8,' changed from ',f20.18,' to ',f20.18)
9060  format (1x,'INFO: total of point ',i8,' is ',f20.18,' not exactly 1.0')

9999  write (*,9000)
      print *, 'while reading '//trim(aname)
      stop 1

end subroutine readit
