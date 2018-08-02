#include "inland_config.h"
#include "inland_compar.h"
!----------------------------------------------------------
subroutine rd_param(irestart)
!----------------------------------------------------------
!  Read various parameters from directory params
!
!---------------------------Code history--------------------------------
! Author: David Price
! 
!----------------------------------------------------------

#ifdef SINGLE_POINT_MODEL
      use inland_comforc ! use everything from comforc.
#endif /* SINGLE_POINT_MODEL */

!      use inland_parameters, only: npoi, npft, nband, ndat, npftu, nsoilay,       &
!                                stemae, rootae, rrootpy, rwoodpy, tempvm_nec,  &
!                                rgrowth, dispuhf, ialoglhf, ialoguhf, avmuir,  &
!                                kfactor, stressfac, indir, datadir
      use inland_control, only: indir, datadir
!      use comage    ! npftu
      use inland_comveg, only: clitrl, clitrm, clitrs, clitll, clitlm, clitls, &
                               clitwl, clitwm, clitws, csoislop, csoislon,     &
                               csoipas, alaiml, alaimu, cgrass, chl, chs, chu, &
                               cleaf, cstem, tblowl, tblows, tblowu, tdripl,   &
                               tdrips, tdripu, wliqlmax, wliqsmax, wliqumax,   &
                               woodnorm, wsnolmax, wsnosmax, wsnoumax, specla, &
                               rhoveg, tauveg, dleaf, dstem, aleaf, aroot,     &
                               awood, optis_params, landusetypemap, landusepftmap, &
                               speclap, awoodp, arootp,aleafp,                 &
                               rhovegvlg, rhovegvlb,      &
                               rhovegvu, rhovegirlg, rhovegirlb, rhovegiru,    &
                               tauvegvlg, tauvegvlb, tauvegvu, tauvegirlg,     &
                               tauvegirlb, tauvegiru

      use inland_comsoi    ! wpudmax
      use inland_parameters
      use inland_compft    ! PFT parameters incl physiological "constants"
      use inland_comtex    ! Soil texture-related parameters
      use inland_combgc    ! Soil biogeochemistry parameters
      use inland_comcrop   ! Crop parameters
      use inland_comnitr   ! N cycling parameters
#ifdef SINGLE_POINT_MODE
      use inland_comveg!,   only: fu,fl
#endif /* SINGLE_POINT_MODEL */
      use inland_comfire, only: blow, bup, betae, Ph, g0, alpha, umax, reparea,&
                                exfire
      use inland_comhour, only:imetyear

      implicit none
  
! Local variables
#ifdef SINGLE_POINT_MODEL
      ! Thiago Veloso (Sep-2012)
      ! Variable below is related to the single-point version
      ! It represents the dimension (number of lines) of the 
      ! input data. It is used later, in subroutine readforc.

!     integer dimforc    ! dimension (number of lines) of input data
#endif /* SINGLE_POINT_MODEL */

      integer*4 parm_unit,& ! file unit assignment for input
                j,i,   & ! PFT number (in range 1 to npftu)
                npft2, & ! number of PFTs reported in params.veg
                npftu2,& ! number of upper canopy PFTs reported in params.veg
                nsoil2,& ! number of soil texture classes reported in params.soi
                !nveg,&  ! number of vegetation classes reported in params.veg - now nvegtype
                optis_i, & ! index in optis_params array
                irestart

      real*8    dummyvar ! use this to read in integers that might be pretending
                         ! to be reals. Also use it as a filler when reading 
                         ! non-existent variables in subroutine ReadItems
      real*8    dummyvars(10)

      character*1024 parm_file
      
      parameter (parm_unit = 9)

#ifdef SINGLE_POINT_MODEL
! Thiago on Set/12: Code below loops over all lines of the input file to
!                   retrieve its dimension, which is saved to variable 
!                   'dimforc'. This procedure avoids both the existence 
!                   of variable 'nbforc' and the declaration of variable 
!                   'dimforc' as a parameter.

       dimforc = 0
       open (unit=20,status='unknown',file=trim(datadir)//'/clim-input')
       do
          read (20,*,end=10)
          dimforc = dimforc + 1
       end do

10     close (20)
       if ( dimforc .eq. 0 ) then
          write (*,*) 'ERROR: could not read clim input file '//trim(datadir)//'/clim-input'
          stop 'make sure INLAND_INDIR is set to proper path'
          stop 1
       end if

#endif /* SINGLE_POINT_MODEL */

      if(imetyear .ne. 9999)then
         nl=0
         open (unit=223,status='unknown',file=trim(datadir)//'/single_agro.txt')
         do
            read (223,*,end=11)
            nl = nl + 1
         end do
11     close (223)

         if ( nl .eq. 0 ) then
             write (*,*) 'ERROR: could not read clim input file '//trim(datadir)//'/single_agro.txt'
             stop 'make sure INLAND_INDIR is set to proper path'
             stop 1
         end if
      endif

! ******************************************************************************
! open the parameter file 'params.can' for input; read in canopy parameters...
     if(isimagro .eq. 0)then
      parm_file = trim(indir)//'/'//'params/canopy'
     else
      parm_file = trim(indir)//'/'//'params/canopy_crop'
     endif
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      call readitem(parm_unit, parm_file, tau15)
      call readitem(parm_unit, parm_file, kc15)
      call readitem(parm_unit, parm_file, ko15)
      call readitem(parm_unit, parm_file, cimax)
      call readitem(parm_unit, parm_file, dummyvar)

      npft2 = nint(dummyvar)
      if (npft2 .ne. npft) then
         write (*, 9003) trim(parm_file), npft2, npft
         goto 9006 ! In the circumstances this seems the best thing to do! 
      end if
9003  format ('rd_param Warning: Number of PFTs in ', A, ' is: ',  &
              I2, ' number in compar.h is: ', I2)  

      call readitem(parm_unit, parm_file, dummyvar)

      npftu2 = nint(dummyvar)
      if (npftu2 .ne. npftu) then
         write (*,9004) trim(parm_file), npftu2, npftu
         goto 9006 ! In the circumstances this seems the best thing to do! 
      end if
9004  format ('rd_param Warning: Number of upper canopy (tree) PFTs ',  &
              'in ', A, ' is: ', I2, ' number in comage.h is: ', I2)  

! Standard photosynthesis parameters for C3 and C4 physiological pathways
      call readitem(parm_unit, parm_file, alpha3)
      call readitem(parm_unit, parm_file, theta3)
      call readitem(parm_unit, parm_file, beta3)
      call readitem(parm_unit, parm_file, betac3)
      call readitem(parm_unit, parm_file, alpha4)
      call readitem(parm_unit, parm_file, theta4)
      call readitem(parm_unit, parm_file, beta4)
      call readitem(parm_unit, parm_file, thetac4)
      call readitem(parm_unit, parm_file, betac4)
      call readitem(parm_unit, parm_file, thetac3)

! Physiological parameters for broadleaf trees
      call readitems(parm_unit, parm_file, 4, gammaub, coefmub, coefbub, &
                     gsubmin, dummyvar, dummyvar, dummyvar, dummyvar,    &
                     dummyvar, dummyvar)

! Physiological parameters for conifer trees
      call readitems(parm_unit, parm_file, 4, gammauc, coefmuc, coefbuc, &
                     gsucmin, dummyvar, dummyvar, dummyvar, dummyvar,    &
                     dummyvar, dummyvar)

! Physiological parameters for shrubs
      call readitems(parm_unit, parm_file, 4, gammals, coefmls, coefbls,       &
                     gslsmin, dummyvar,dummyvar, dummyvar, dummyvar, dummyvar, &
                     dummyvar)

! Physiological parameters for C4 grasses
      call readitems(parm_unit, parm_file, 4, gammal4, coefml4, coefbl4,       &
                     gsl4min, dummyvar,dummyvar, dummyvar, dummyvar, dummyvar, &
                     dummyvar)

! Physiological parameters for C3 grasses
      call readitems(parm_unit, parm_file, 4, gammal3, coefml3, coefbl3,       &
                     gsl3min, dummyvar,dummyvar, dummyvar, dummyvar, dummyvar, &
                     dummyvar)

! Physiological parameters for C3 crops
      call readitems(parm_unit, parm_file, 4,                   &
               gammac3, coefmc3, coefbc3, gsc3min, dummyvar,    &
               dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

! Physiological parameters for C4 crops
      call readitems(parm_unit, parm_file, 4,                   &
               gammac4, coefmc4, coefbc4, gsc4min, dummyvar,    &
               dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

!!! DTP 2001/06/05: Note that I have included a new variable here: vmax_pft
!                   This reads in values of vmax assigned initially to each
!                   PFT. These values are then transferred to vmaxub, vmaxuc,
!                   vmaxls, vmaxl4, vmaxl3, in the modified PHYSIOLOGY.F 
!                   module (subroutine STOMATA).

! Castanho HP, 2013 included maps vegetation properties (vmax,specla,awood,aroot,aleaf,tauwood) in input file read in io.f,
! we are keeping the constant homogeneous parameters in params.can and changing the names adding 'p' 
   if(isimagro.gt.0) then
    if (irestart .eq. 0)then
      do j = 1, npft 
         call readitems(parm_unit, parm_file, 8, vmax_pftp(j), speclap(j), &
                        tauleaf(j), tauroot(j), tauwood0p(j),aleafp(j),    &    	!Castanho HP, 2013 commented aleafp(j) 
                        arootp(j), awoodp(j), dummyvar, dummyvar, dummyvar) 
      end do
    else
      do j = 1, npft 
         call readitems(parm_unit, parm_file, 4, vmax_pftp(j), speclap(j), &
                        tauleaf(j), tauroot(j), tauwood0p(j), dummyvar,    &
                        dummyvar, dummyvar) 
      end do
     endif
    else
      do j = 1, npft 
            call readitems(parm_unit, parm_file, 8, vmax_pftp(j), speclap(j), &
                        tauleaf(j), tauroot(j), tauwood0p(j),aleafp(j),    &   
                        arootp(j), awoodp(j), dummyvar, dummyvar, dummyvar) 
      end do
    endif

      call readitem(parm_unit, parm_file, woodnorm)

      if(isimagro.gt.0) then

         call readitem(parm_unit, parm_file, rhovegvlg)
         call readitem(parm_unit, parm_file, rhovegvlb)
         call readitem(parm_unit, parm_file, rhovegvu)
         call readitem(parm_unit, parm_file, rhovegirlg)
         call readitem(parm_unit, parm_file, rhovegirlb)
         call readitem(parm_unit, parm_file, rhovegiru)
         call readitem(parm_unit, parm_file, tauvegvlg)
         call readitem(parm_unit, parm_file, tauvegvlb)
         call readitem(parm_unit, parm_file, tauvegvu)
         call readitem(parm_unit, parm_file, tauvegirlg)
         call readitem(parm_unit, parm_file, tauvegirlb)
         call readitem(parm_unit, parm_file, tauvegiru)
      
      endif

      do j = 1, nband
         call readitems(parm_unit, parm_file, 2, rhoveg(j,1), rhoveg(j,2), &
                        dummyvar, dummyvar, dummyvar, dummyvar, dummyvar,  &
                        dummyvar, dummyvar, dummyvar) 
      end do

      do j = 1, nband
         call readitems(parm_unit, parm_file, 2, tauveg(j,1), tauveg(j,2), &
                        dummyvar, dummyvar, dummyvar, dummyvar, dummyvar,  &
                        dummyvar, dummyvar, dummyvar)   
      end do
      
      if(isimagro .gt. 0)then
!         call readitem(parm_unit, parm_file, chifuz)
!         call readitem(parm_unit, parm_file, chiflz)
      else 
         call readitem(parm_unit, parm_file, chifuz)
         call readitem(parm_unit, parm_file, chiflz)
      endif

      call readitems(parm_unit, parm_file, 2, dleaf(1), dleaf(2), dummyvar,    &
                     dummyvar, dummyvar, dummyvar,dummyvar, dummyvar, dummyvar,&
                     dummyvar)

      call readitems(parm_unit, parm_file, 2, dstem(1), dstem(2), dummyvar,    &
                     dummyvar, dummyvar, dummyvar,dummyvar, dummyvar, dummyvar,&
                     dummyvar)

      call readitem(parm_unit, parm_file, alaimu)
      call readitem(parm_unit, parm_file, alaiml)

      call readitem(parm_unit, parm_file, cleaf)
      call readitem(parm_unit, parm_file, cstem)
      call readitem(parm_unit, parm_file, cgrass)

      call readitem(parm_unit, parm_file, chs)
      call readitem(parm_unit, parm_file, chu)
      call readitem(parm_unit, parm_file, chl)

      call readitem(parm_unit, parm_file, wliqumax)
      call readitem(parm_unit, parm_file, wliqsmax)
      call readitem(parm_unit, parm_file, wliqlmax)

      call readitem(parm_unit, parm_file, wsnoumax)
      call readitem(parm_unit, parm_file, wsnosmax)
      call readitem(parm_unit, parm_file, wsnolmax)

      call readitem(parm_unit, parm_file, tdripu)
      call readitem(parm_unit, parm_file, tdrips)
      call readitem(parm_unit, parm_file, tdripl)

      call readitem(parm_unit, parm_file, tblowu)
      call readitem(parm_unit, parm_file, tblows)
      call readitem(parm_unit, parm_file, tblowl)

! Canopy parameters identified on Optis' efforts - fzm
      call readitem(parm_unit, parm_file, stemae)
      call readitem(parm_unit, parm_file, rootae)
      call readitem(parm_unit, parm_file, rrootpy)
      call readitem(parm_unit, parm_file, rwoodpy)
      call readitem(parm_unit, parm_file, tempvm_nec)
      call readitem(parm_unit, parm_file, rgrowth)
      call readitem(parm_unit, parm_file, dispuhf)
      call readitem(parm_unit, parm_file, ialoglhf)
      call readitem(parm_unit, parm_file, ialoguhf)
      call readitem(parm_unit, parm_file, avmuir)
#ifdef SINGLE_POINT_MODE
      call readitem(parm_unit, parm_file, fu)
      call readitem(parm_unit, parm_file, fl)
#endif /* SINGLE_POINT_MODEL */
      close (parm_unit)

! ******************************************************************************
! open the parameter file 'params.veg' for input; read in vegetation PFT
! parameters...
     if(isimagro .eq. 0)then
      parm_file = trim(indir)//'/'//'params/vegetation'
     else
      parm_file = trim(indir)//'/'//'params/vegetation_crop'
     endif
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      do j = 1, npft
         call readitems(parm_unit, parm_file, 4, TminL(j), TminU(j), Twarm(j), &
                        GDD(j),dummyvar, dummyvar, dummyvar, dummyvar,         &
                        dummyvar, dummyvar)
      end do 

      call readitem(parm_unit, parm_file, dummyvar)
      nvegtype = nint(dummyvar)

!   There's no better place to allocate plai_init than here, as here we know its
! dimension size and we use it right away.
      allocate(plai_init(4,nvegtype))
      plai_init(:,:) = 0. ! do we really need to zero it here?

      do j = 1, nvegtype      
         call readitems(parm_unit, parm_file, 4, plai_init(1,j),          &
                        plai_init(2,j), plai_init(3,j), plai_init(4,j),   &
                        dummyvar, dummyvar, dummyvar, dummyvar, dummyvar, &
                        dummyvar)
      end do

      call readitem(parm_unit, parm_file, plaiupper)
      call readitem(parm_unit, parm_file, plailower)
      call readitem(parm_unit, parm_file, xminlai)
      call readitem(parm_unit, parm_file, sapfrac_init)
      call readitem(parm_unit, parm_file, beta1)
      call readitem(parm_unit, parm_file, beta2)

! Parameters identified in Optis' efforts - fzm and et
! Read the values into optis_params, and propagate them in inland_alloc, 
! because for now we don't know npoi
      do optis_i = 1, 12
         call readitem(parm_unit, parm_file, optis_params(optis_i))
      end do
! TODO - what happened to csoislo?
#ifndef SINGLE_POINT_MODEL
! new parameters for anthropic vegtypes

      call readitem(parm_unit, parm_file, dummyvar)
      nlandusetype = nint(dummyvar)

      allocate(landusetypemap(nvegtype))
      landusetypemap(:) = 0
      do j = 1, nvegtype
         call readitem(parm_unit, parm_file, dummyvar)
         landusetypemap(j) = nint(dummyvar)
      end do

      allocate(landusepftmap(nlandusetype,npft))
      landusepftmap(:,:) = 0.
      do j = 1, npft
         call readitems(parm_unit, parm_file, 4,dummyvars(1),dummyvars(2), &
                        dummyvars(3),dummyvars(4),dummyvars(5),dummyvars(6), &
                        dummyvars(7),dummyvars(8),dummyvars(9),dummyvars(10))
         do i = 1, nlandusetype
            landusepftmap(i,j) = nint(dummyvars(i))
         end do
      end do

#endif  /* SINGLE_POINT_MODEL */
      close (parm_unit)
      

! ******************************************************************************
! open the parameter file 'params.soi' for input; read in soil parameters...
     if(isimagro .eq. 0)then
      parm_file = trim(indir)//'/'//'params/soil'
     else
      parm_file = trim(indir)//'/'//'params/soil_crop'
     endif
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      call readitem(parm_unit, parm_file, dummyvar)
      nsoilay = nint(dummyvar)

      call inland_inneralloc

      do j = 1, nsoilay
         call readitem(parm_unit, parm_file, hsoi(j))
      end do

      call readitem(parm_unit, parm_file, dummyvar)
      nslaym = nint(dummyvar)
      call readitem(parm_unit, parm_file, bperm)
      call readitem(parm_unit, parm_file, wpudmax)
      call readitem(parm_unit, parm_file, zwpmax)

      call readitem(parm_unit, parm_file, dummyvar)

      nsoil2 = nint(dummyvar)

      if (nsoil2 .ne. ndat) then ! FIXME: Is this needed????
         write (*, 9031) trim(parm_file), nsoil2
         write (*, 9032) ndat
         goto 9006 ! In the circumstances this seems the best thing to do! 
      end if

9031  format (' rd_param Warning: Number of soil types in ', A10, ' is: ', I2)
9032  format (' Number of soil types in comtex.h is: ', I2)  

      do j = 1, ndat
         call readitems(parm_unit, parm_file, 10, texdat(1,j), texdat(2,j),   &
                        texdat(3,j), porosdat(j), sfielddat(j), swiltdat(j), &
                        bexdat(j), suctiondat(j), hydrauldat(j), cpwfdat(j))      ! Kai included cpwfdat
      end do

      call readitem(parm_unit, parm_file, lig_frac)
      call readitem(parm_unit, parm_file, fbsom)
      call readitem(parm_unit, parm_file, effac)

      call readitems(parm_unit, parm_file, 8, cnr(1), cnr(2), cnr(3), cnr(4), &
                     cnr(5), cnr(6), cnr(7), cnr(8), dummyvar, dummyvar)

    if (irestart .eq. 0)then
      call readitems(parm_unit, parm_file, 6, fmax, rconst, cnleaf,cnwood , &
                     cnroot, h20, dummyvar, dummyvar, dummyvar, &
                     dummyvar,dummyvar)
      call readitems(parm_unit, parm_file, 5, co, cf, fnitrate, nloss, cndepth, &
                     dummyvar, dummyvar, dummyvar,dummyvar,dummyvar)
    else
      call readitems(parm_unit, parm_file, 2, fmax, rconst, &
                     dummyvar, dummyvar, dummyvar, &
                     dummyvar,dummyvar)
      call readitems(parm_unit, parm_file, 4, co, cf, fnitrate, nloss, &
                     dummyvar, dummyvar, dummyvar,dummyvar,dummyvar)
    endif

      call readitems(parm_unit, parm_file, 9, klm, kls, kll, krm, krs, krl, &
                     kwm, kws, kwl, dummyvar)

      call readitems(parm_unit, parm_file, 7, kbn, kbp, knb, kns, kpb, kps, &
                     ksb, dummyvar, dummyvar, dummyvar)

      call readitems(parm_unit, parm_file, 9, ylm, yrm, ywm, yls, yrs, yws, &
                     yll, yrl, ywl, dummyvar)

      call readitems(parm_unit, parm_file, 7, ybn, ybp, yps, yns, ysb, ypb, &
                     ynb, dummyvar, dummyvar, dummyvar)

! Parameters identified in Optis' efforts
      call readitem(parm_unit, parm_file, stressfac)

! We read wsoi value and add it to optis_params
      call readitem(parm_unit, parm_file, optis_params(13))

      call readitem(parm_unit, parm_file, kfactor)

      close (parm_unit)
! ******************************************************************************
! open the parameter file 'params.fir' for input; read in fire parameters...
      parm_file = trim(indir)//'/'//'params/fire'
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001) 

      call readitem(parm_unit, parm_file, blow)
      call readitem(parm_unit, parm_file, bup)
      call readitem(parm_unit, parm_file, betae)
      call readitem(parm_unit, parm_file, Ph)
      call readitem(parm_unit, parm_file, g0)
      call readitem(parm_unit, parm_file, alpha)
      call readitem(parm_unit, parm_file, umax)
      call readitem(parm_unit, parm_file, reparea)
      call readitem(parm_unit, parm_file, exfire)

! ******************************************************************************

#ifndef SINGLE_POINT_MODEL
! ******************************************************************************
! open the parameter file 'params.crp' for input; read in crop parameters...
      parm_file = trim(indir)//'/'//'params/crop'
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      do j = scpft, ecpft
         call readitems(parm_unit, parm_file, 8, lotemp(j), hitemp(j), &
                        drought(j), f1(j), f2(j), baset(j), mxtmp(j),  &
                        tkill(j), dummyvar, dummyvar) 
      end do

      do j = scpft, ecpft 
         call readitems(parm_unit, parm_file, 10, laicons(j), allconsl(j),     &
                        allconss(j), laimx(j), arooti(j), arootf(j), aleaff(j), &
                        astemf(j), declfact(j), fleafi(j)) 
      end do

      do j = scpft, ecpft 
         call readitems(parm_unit, parm_file, 6, ptemp(j), pmintemp(j), &
                        pmmin(j), pdmin(j), pcm(j), pcd(j), dummyvar,   &
                        dummyvar, dummyvar, dummyvar)
      end do

      do j = scpft, ecpft 
         call readitems(parm_unit, parm_file, 8, hybgdd(j), gddmin(j), &
                        mxgddgf(j), mxdgfi(j), mxmat(j), lfemerg(j),   &
                        grnfill(j), bfact(j), dummyvar, dummyvar) 
      end do

        pcm(16)=mod(pmmin(16)+(mxmat(16)/30),12.)+1
        pcd(16)=pdmin(16)
      
      do j = 1, 2 ! number of wheat crop types 
         call readitems(parm_unit, parm_file, 7, grnwht(j), fleafiwht(j), &
                        mgddgf(j), mxmatwht(j), fnlfmxw(j), fngrmxw(j),   &
                        fnoptw(j), dummyvar, dummyvar, dummyvar) 
      end do

      call readitem(parm_unit, parm_file, ztopmxsoy)
      call readitem(parm_unit, parm_file, ztopmxwht)
      call readitem(parm_unit, parm_file, ztopmxmze)
      call readitem(parm_unit, parm_file, ztopmxsgc)
      call readitem(parm_unit, parm_file, dummyvar)
      nratoon=int(dummyvar)
      call readitem(parm_unit, parm_file, alphac)
      call readitem(parm_unit, parm_file, gnmin)
      call readitem(parm_unit, parm_file, smax)
      call readitem(parm_unit, parm_file, availn)
      call readitem(parm_unit, parm_file, cnmax)
      call readitem(parm_unit, parm_file, rootd)
      call readitem(parm_unit, parm_file, laidcs)
      call readitem(parm_unit, parm_file, mxtsuc)
      call readitem(parm_unit, parm_file, sf1)
      call readitem(parm_unit, parm_file, ipf1)
      call readitem(parm_unit, parm_file, ecf2)
      call readitem(parm_unit, parm_file, ipf2)
      call readitem(parm_unit, parm_file, sf3)
      call readitem(parm_unit, parm_file, ipf3)
      call readitem(parm_unit, parm_file, ecf4)
      call readitem(parm_unit, parm_file, ipf4)
      call readitem(parm_unit, parm_file, ecf5)
      call readitem(parm_unit, parm_file, tf5)
      call readitem(parm_unit, parm_file, wf5)
      call readitem(parm_unit, parm_file, ecf6)
      call readitem(parm_unit, parm_file, ipf6)
      call readitem(parm_unit, parm_file, laidc)
      call readitem(parm_unit, parm_file, ldf)
      call readitem(parm_unit, parm_file, tmld)
      call readitem(parm_unit, parm_file, ecf7)
      call readitem(parm_unit, parm_file, firecane)
     
      do j = scpft, ecpft 
         call readitems(parm_unit, parm_file, 10, cgrain(j), convfact(j),    &
                        maxhi(j), fyield(j), cfrac(j), fnlfmx(j), fngrmx(j), &
                        sratio(j), rratio(j), fnopt(j)) 
      end do

      do j = scpft, ecpft 
         call readitems(parm_unit, parm_file, 1, grainmoisture(j)) 
      end do

      close (parm_unit)
#endif  /* SINGLE_POINT_MODEL */
! ******************************************************************************
! open the parameter file 'params.hyd' for input; read in vegetation PFT parameters...
!      parm_file = 'params.hyd' ! What is this file supposed to contain???
!      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)
!      close (parm_unit) 
! ******************************************************************************

      write (*,*)'rd_param: All data read in from parameter files successfully.'
      return ! subroutine rd_param

! Error handling targets
9001  write (*,2001) TRIM(parm_file), parm_unit 
2001  format ('rd_param: Error opening parameter file "', A,'" on unit ', I2)      
      stop

9006  write (*,2006) TRIM(parm_file), parm_unit
2006  format ('rd_param: Data inconsistency in parameter file "',A, &
              '" on unit ', I2) 
      close (parm_unit)
      stop
end subroutine rd_param
!*******************************************************************************
