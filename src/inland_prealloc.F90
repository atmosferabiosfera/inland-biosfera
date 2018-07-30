#include "inland_config.h"
! Inland's Pre allocate phase subroutine

!   This subroutine Allocates variables as they are needed by rd_param
! and readit subroutines.
subroutine inland_prealloc
      use inland_parameters, only: npft,nband,nsoilay,ndat,plona,plata,nlonsub, &
                                 nlatsub
      use inland_compft, only: vmax_pft,tauleaf,tauroot,tauwood0,TminL,TminU, &
                               Twarm,GDD,lotemp,hitemp,drought,f1,f2, tauwood0p,vmax_pftp
!ver depois
      use inland_comveg, only: specla,aleaf,aroot,awood,rhoveg,tauveg,dleaf,  &
                               dstem, clitrl, clitrm, clitrs, clitll, clitlm, &
                               clitls, clitwl, clitwm, clitws, csoislon,      &
                               csoislop, csoipas, fu, fl, speclap, awoodp,    &
                               arootp, aleafp
      use inland_comtex ! only: all!
      use inland_combgc, only: cnr
      use inland_comsoi, only: hsoi !, sand, clay, wsoi
      use inland_comwork ! only: all!
      use inland_combcs, only: xintopo,xinveg,deltat,xint,clmwet, &
                               clmtrng,clmt,clmprec,clmcld,clmq,clmwind, &
                               xintrng,xinprec,xincld,xinq,xinwet,xinwetmon, &
                               xinprecmon, xintrngmon, xincldmon, xinqmon, &
                               xinwindmon, xintmon,obswet,obstrng,obst,obsprec, &
                               obscld, obsq, obswind,xinwind
      use inland_comcrop

      implicit none

      integer*1 bytevar ! variables for testing ocean (netcdf nodata) values
      integer*2 shortvar
      integer   intvar
      real*8    realvar
      real*4    real4var

! Castanho HP, 2013 included dimensions npoi (i) in aleaf, awood, aroot, tauwood, specla, vmax when appropriate 

!! Variables for rd_param
!==========================
! Allocate variables in comcrop
      allocate(baset(npft),mxtmp(npft),tkill(npft),laicons(npft),allconsl(npft), &
               allconss(npft),laimx(npft),arooti(npft),arootf(npft),aleaff(npft), &
               astemf(npft),declfact(npft),fleafi(npft),hybgdd(npft),gddmin(npft),&                     
               lfemerg(npft),grnfill(npft),mxgddgf(npft),mxdgfi(npft),mxmat(npft),&
               bfact(npft),arepr(npft),astemi(npft),aleafi(npft),fleaf(npft))

      baset(:) = 0.
      mxtmp(:) = 0.
      tkill(:) = 0.
      laicons(:) = 0.
      allconsl(:) = 0.
      allconss(:) = 0.
      laimx(:) = 0.
      arooti(:) = 0.
      arootf(:) = 0.
      aleaff(:) = 0.
      astemf(:) = 0.
      declfact(:) = 0.
      fleafi(:) = 0.
      hybgdd(:) = 0.
      gddmin(:) = 0.   
      lfemerg(:) = 0.  
      grnfill(:) = 0.  
      mxgddgf(:) = 0.  
      mxdgfi(:) = 0.   
      mxmat(:) = 0.	 
      bfact(:) = 0.
      arepr(:) = 0.
      astemi(:) = 0.   
      aleafi(:) = 0.   
      fleaf(:) = 0.

      allocate(ptemp(npft),pmintemp(npft),pmmin(npft),pdmin(npft),pcm(npft),pcd(npft))

      ptemp(:) = 0.
      pmintemp(:) = 0.
      pmmin(:) = 0.
      pdmin(:) = 0.
      pcm(:) = 0.
      pcd(:) = 0.

      allocate(grnwht(2),fleafiwht(2),mgddgf(2),mxmatwht(2),fnlfmxw(2),fngrmxw(2),&
               fnoptw(2))

      grnwht(:) = 0.
      fleafiwht(:) = 0.
      mgddgf(2) = 0.
      mxmatwht(2) = 0.
      fnlfmxw(2) = 0.
      fngrmxw(2) = 0.
      fnoptw(2) = 0.

      allocate (cgrain(npft),convfact(npft),maxhi(npft),fyield(npft),cfrac(npft),&
                fnlfmx(npft),fngrmx(npft),sratio(npft),rratio(npft),fnopt(npft),grainmoisture(npft))

      cgrain(:) = 0.
      convfact(:) = 0.
      maxhi(:) = 0.
      fyield(:) = 0.
      cfrac(:) = 0.
      fnlfmx(:) = 0.
      fngrmx(:) = 0.
      sratio(:) = 0.
      rratio(:) = 0.
      fnopt(:) = 0.
      grainmoisture(:) = 0.
	     
! From compft

! Castanho HP, 2013 moved vmax_pft and tauwood0 2 dimension to allocate in readit 
! Castanho HP, 2013 added *p to one dimension that is read in rd param

      allocate(tauleaf(npft),tauroot(npft), tauwood0p(npft), vmax_pftp(npft), &  
               TminL(npft),TminU(npft),Twarm(npft),GDD(npft),lotemp(npft),&
               hitemp(npft),drought(npft),f1(npft),f2(npft))

      tauleaf(:) = 0.
      tauroot(:) = 0.
!     tauwood0(:,:) = 0.			! Castanho HP, 2013
      tauwood0p(:) = 0.				! Castanho HP, 2013
!     vmax_pft(:,:) = 0.			! Castanho HP, 2013
      vmax_pftp(:) = 0.				! Castanho HP, 2013
      TminL(:) = 0.
      TminU(:) = 0.
      Twarm(:) = 0.
      GDD(:) = 0.
      lotemp(:) = 0.
      hitemp(:) = 0.
      drought(:) = 0.
      f1(:) = 0.
      f2(:) = 0.

! From comveg

! Castanho HP, 2013 keep one dimension allocation addded p names

      allocate(speclap(npft),arootp(npft),awoodp(npft), aleafp(npft),  &
               rhoveg(nband,2),tauveg(nband,2))

      speclap(:) = 0.				! Castanho HP, 2013
      aleafp(:)  = 0.				! Castanho HP, 2013 function of the others
      arootp(:)  = 0.				! Castanho HP, 2013
      awoodp(:)  = 0.				! Castanho HP, 2013
      rhoveg(:,:) = 0.
      tauveg(:,:) = 0.

      allocate(dleaf(2),dstem(2))

      dleaf(:) = 0.
      dstem(:) = 0.

! Allocate variables in comsoi
      allocate(hsoi(nsoilay+1))! FIXME: really needed to allocate nsoilay+1?
      hsoi(:) = 0.

! Allocate variables in comtex
      allocate(texdat(3,ndat),porosdat(ndat),sfielddat(ndat),swiltdat(ndat), &
               bexdat(ndat),suctiondat(ndat),hydrauldat(ndat),cpwfdat(ndat))
      texdat(:,:) = 0.     
      porosdat(:) = 0.
      sfielddat(:) = 0.
      swiltdat(:) = 0.
      bexdat(:) = 0.
      suctiondat(:) = 0.
      hydrauldat(:) = 0.
      cpwfdat(:) = 0.                         ! Kai included cpwfdat

! Allocate variables in combgc
      allocate(cnr(10))
      cnr(:) = 0.

!! Variables for readit
!================================
! Allocate variables in comwork
!--------------------------------

! Get comwork's ndim4 size
!    The base values are now in inland_parameters.F90 module; 
!    in the future they will be read from inland.infile
      ndim2 = nlon*nlat
      ndim4 = max(nlon,nlat,nband,nsoilay,nsnolay,npft)
      ndim3=nlon*nlat*max(nband,nsoilay,nsnolay,npft)
      if ( mlpt .gt. 1 ) ndim3 = ndim3*(mlpt+1)

      allocate(lonscale(nlon),latscale(nlat))
      allocate(work(ndim2))
      allocate(cdummy(ndim3),cdummyint(ndim3))
      lonscale(:) = 0.
      latscale(:) = 0.
      work(:) = 0.
      cdummy(:) = 0.
      cdummyint(:) = 0

! removed allocation of variables with npoi dimension, now in readit

end subroutine inland_prealloc
