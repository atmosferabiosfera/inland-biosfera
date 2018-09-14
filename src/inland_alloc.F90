#include "inland_config.h"
!TODO: header, make allocation subroutines inside each corresponding module
! ---------------------------------------------------------------------
subroutine alloc(irestart)
!-----------------------------------------------------------------------
!     Allocate variables
!---------------------------Code history--------------------------------
!     Original version: A.A. Mirin, March 2001
! ---------------------------------------------------------------------
      use inland_com1d
      use inland_lsmmapib
      use inland_parameters
      use inland_control
      use inland_comatm
      use inland_combcs
      use inland_comfire
      use inland_comcrop
      use inland_comhour
      use inland_comhyd
      use inland_comnitr
      use inland_compft
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comwork, only: seedvec
      use inland_combgc, only:cnroot
#ifdef SINGLE_POINT_MODEL
      use inland_comforc ! use everything from comforc.
#endif

      implicit none

      integer irestart

! Variables in comveg identified on optis' parameter set, was in inland_prealloc
      allocate(clitrl(lbeg:lend),clitrm(lbeg:lend),clitrs(lbeg:lend), &
               clitll(lbeg:lend),clitlm(lbeg:lend),clitls(lbeg:lend), &
               clitwl(lbeg:lend),clitwm(lbeg:lend),clitws(lbeg:lend), &
               csoislon(lbeg:lend),csoislop(lbeg:lend),csoipas(lbeg:lend))
      clitrl(:) = optis_params(1)
      clitrm(:) = optis_params(2)
      clitrs(:) = optis_params(3)
      clitll(:) = optis_params(4)
      clitlm(:) = optis_params(5)
      clitls(:) = optis_params(6)
      clitwl(:) = optis_params(7)
      clitwm(:) = optis_params(8)
      clitws(:) = optis_params(9)
! TODO - what happened to csoislo?
      csoislop(:) = optis_params(10)
      csoislon(:) = optis_params(11)
      csoipas(:) = optis_params(12)

      allocate(wsoi(lbeg:lend,nsoilay))
      wsoi(:,:) = optis_params(13)

! Zero constants in inland_control part of inland_parameters (subroutine)
      stef = 0.
      vonk = 0.
      grav = 0.
      tmelt = 0.
      hfus = 0.
      hvap = 0.
      hsub = 0.
      ch2o = 0.
      cice = 0.
      cair = 0.
      cvap = 0.
      rair = 0.
      rvap = 0.
      cappa = 0.
      rhow = 0.

! Allocate variables in inland_control part of inland_parameters (subroutine)
! Please note, 'mpt' is the equivalent of 'npoi' for a single MPI process.
! FIXME: Decide between mpt or lbeg:lend!
      allocate(vzero(lbeg:lend))

! Zero (non-namelist) constants in inland_control.F90
      spincons = 0.
      eqyears = 0
      spinfrac = 0.
      nspinsoil = 0
      iyear = 0
      imonth = 0
      if(irestart .eq. 0)then
         iday = 0
      endif

! Allocate variables in lsmmapib
!     (begpnt and numpnt allocated separately)
      allocate(ixy(lbeg:lend),jxy(lbeg:lend),lati(lbeg:lend),loni(lbeg:lend), &
               procj(lbeg:lend),procn(lbeg:lend))

      ixy(:) = 0
      jxy(:) = 0
      lati(:) = 0.
      loni(:) = 0.
      procj(:) = 0
      procn(:) = 0

      allocate(ixy1(npoi1),jxy1(npoi1),lati1(npoi1),loni1(npoi1),procj1(npoi1),&
               procn1(npoi1))

      ixy1(:) = 0
      jxy1(:) = 0
      lati1(:) = 0.
      loni1(:) = 0.
      procj1(:) = 0
      procn1(:) = 0

! Allocate variables in com1d
      allocate(terml(lbeg:lend,7),termu(lbeg:lend,7))
      allocate(scalcoefl(lbeg:lend,4),scalcoefu(lbeg:lend,4))
      allocate(firb (lbeg:lend),&
               firg (lbeg:lend),&
               firi (lbeg:lend),&
               firl (lbeg:lend),&
               firs (lbeg:lend),&
               firu (lbeg:lend),&
               fsena (lbeg:lend),&
               fseng (lbeg:lend),&
               fseni(lbeg:lend),&
               fsenl (lbeg:lend),&
               fsens (lbeg:lend),&
               fsenu (lbeg:lend),&
               fvapa (lbeg:lend),&
               fvapg (lbeg:lend),&
               fvapi (lbeg:lend),&
               fvaplt (lbeg:lend),&
               fvaplw (lbeg:lend),&
               fvaps (lbeg:lend),&
               fvaput (lbeg:lend),&
               fvapuw (lbeg:lend),&
               raing(lbeg:lend),&
               rainl(lbeg:lend),&
               rainu(lbeg:lend),&
               snowg(lbeg:lend),&
               snowl(lbeg:lend),&
               snowu(lbeg:lend),&
               soli (lbeg:lend),&
               solg (lbeg:lend),&
               soll (lbeg:lend),&
               sols (lbeg:lend),&
               solu (lbeg:lend),&
               traing(lbeg:lend),&
               trainl(lbeg:lend),&
               trainu(lbeg:lend),&
               tsnowg(lbeg:lend),&
               tsnowl(lbeg:lend),&
               tsnowu(lbeg:lend))
      allocate(nsol(numlv))
      allocate(abupd(lbeg:lend),&
               abupi(lbeg:lend),&
               ablod(lbeg:lend),&
               abloi(lbeg:lend),&
               albsnd(lbeg:lend),&
               albsni(lbeg:lend),&
               albsod(lbeg:lend),&
               albsoi(lbeg:lend),&
               dummy(lbeg:lend),&
               flodd(lbeg:lend),&
               flodi(lbeg:lend),&
               floii(lbeg:lend),&
               fupdd(lbeg:lend),&
               fupdi(lbeg:lend),&
               fupii(lbeg:lend),&
               relod(lbeg:lend),&
               reloi(lbeg:lend),&
               reupd(lbeg:lend),&
               reupi(lbeg:lend),&
               sol2d(lbeg:lend),&
               sol2i(lbeg:lend),&
               sol3d(lbeg:lend),&
               sol3i(lbeg:lend))
      allocate(indsol(numlv,lend-lbeg+1))
      allocate(aloga(lbeg:lend),&
               alogav(lbeg:lend),&
               alogg(lbeg:lend),&
               alogi(lbeg:lend),&
               alogl(lbeg:lend),&
               alogu(lbeg:lend),&
               alog1(lbeg:lend),&
               alog2(lbeg:lend),&
               alog3(lbeg:lend),&
               alog4(lbeg:lend),&
               bdl(lbeg:lend),&
               bdu(lbeg:lend),&
               cl(lbeg:lend),&
               cp(lbeg:lend),&
               cu(lbeg:lend),&
               dil(lbeg:lend),&
               displ(lbeg:lend),&
               dispu(lbeg:lend),&
               diu(lbeg:lend),&
               exphl(lbeg:lend),&
               exphu(lbeg:lend),&
               expl(lbeg:lend),&
               expu(lbeg:lend),&
               pfluxl(lbeg:lend),&
               pfluxs(lbeg:lend),&
               pfluxu(lbeg:lend),&
               rhoa(lbeg:lend),&
               richl(lbeg:lend),&
               richu(lbeg:lend),&
               sg(lbeg:lend),&
               si(lbeg:lend),&
               strahl(lbeg:lend),&
               strahu(lbeg:lend),&
               straml(lbeg:lend),&
               stramu(lbeg:lend),&
               u1(lbeg:lend),&
               u12(lbeg:lend),&
               u2(lbeg:lend),&
               u3(lbeg:lend),&
               u34(lbeg:lend),&
               u4(lbeg:lend),&
               za(lbeg:lend),&
               z1(lbeg:lend),&
               z12(lbeg:lend),&
               z2(lbeg:lend),&
               z3(lbeg:lend),&
               z34(lbeg:lend),&
               z4(lbeg:lend),&
               taux(lbeg:lend),&
               tauy(lbeg:lend),&
               ts2(lbeg:lend),&
               qs2(lbeg:lend),&
               tfaca(lbeg:lend))
      allocate(fwetl (lbeg:lend),&
               fwets (lbeg:lend),&
               fwetu (lbeg:lend),&
               fwetlx (lbeg:lend),&
               fwetsx (lbeg:lend),&
               fwetux (lbeg:lend),&
               rliql (lbeg:lend),&
               rliqs (lbeg:lend),&
               rliqu (lbeg:lend))

      terml(:,:) = 0.
      termu(:,:) = 0.
      scalcoefl(:,:) = 0.
      scalcoefu(:,:) = 0.
      firb(:) = 0.
      firg(:) = 0.
      firi(:) = 0.
      firl(:) = 0.
      firs(:) = 0.
      firu(:) = 0.
      fsena(:) = 0.
      fseng(:) = 0.
      fseni(:) = 0.
      fsenl(:) = 0.
      fsens(:) = 0.
      fsenu(:) = 0.
      fvapa(:) = 0.
      fvapg(:) = 0.
      fvapi(:) = 0.
      fvaplt(:) = 0.
      fvaplw(:) = 0.
      fvaps(:) = 0.
      fvaput(:) = 0.
      fvapuw(:) = 0.
      raing(:) = 0.
      rainl(:) = 0.
      rainu(:) = 0.
      snowg(:) = 0.
      snowl(:) = 0.
      snowu(:) = 0.
      soli(:) = 0.
      solg(:) = 0.
      soll(:) = 0.
      sols(:) = 0.
      solu(:) = 0.
      traing(:) = 0.
      trainl(:) = 0.
      trainu(:) = 0.
      tsnowg(:) = 0.
      tsnowl(:) = 0.
      tsnowu(:) = 0.
      nsol(:) = 0
      abupd(:) = 0.
      abupi(:) = 0.
      ablod(:) = 0.
      abloi(:) = 0.
      albsnd(:) = 0.
      albsni(:) = 0.
      albsod(:) = 0.
      albsoi(:) = 0.
      dummy(:) = 0.
      flodd(:) = 0.
      flodi(:) = 0.
      floii(:) = 0.
      fupdd(:) = 0.
      fupdi(:) = 0.
      fupii(:) = 0.
      relod(:) = 0.
      reloi(:) = 0.
      reupd(:) = 0.
      reupi(:) = 0.
      sol2d(:) = 0.
      sol2i(:) = 0.
      sol3d(:) = 0.
      sol3i(:) = 0.
      indsol(:,:) = 0
      aloga(:) = 0.
      alogav(:) = 0.
      alogg(:) = 0.
      alogi(:) = 0.
      alogl(:) = 0.
      alogu(:) = 0.
      alog1(:) = 0.
      alog2(:) = 0.
      alog3(:) = 0.
      alog4(:) = 0.
      bdl(:) = 0.
      bdu(:) = 0.
      cl(:) = 0.
      cp(:) = 0.
      cu(:) = 0.
      dil(:) = 0.
      displ(:) = 0.
      dispu(:) = 0.
      diu(:) = 0.
      exphl(:) = 0.
      exphu(:) = 0.
      expl(:) = 0.
      expu(:) = 0.
      pfluxl(:) = 0.
      pfluxs(:) = 0.
      pfluxu(:) = 0.
      rhoa(:) = 0.
      richl(:) = 0.
      richu(:) = 0.
      sg(:) = 0.
      si(:) = 0.
      strahl(:) = 0.
      strahu(:) = 0.
      straml(:) = 0.
      stramu(:) = 0.
      u1(:) = 0.
      u12(:) = 0.
      u2(:) = 0.
      u3(:) = 0.
      u34(:) = 0.
      u4(:) = 0.
      za(:) = 0.
      z1(:) = 0.
      z12(:) = 0.
      z2(:) = 0.
      z3(:) = 0.
      z34(:) = 0.
      z4(:) = 0.
      taux(:) = 0.
      tauy(:) = 0.
      ts2(:) = 0.
      qs2(:) = 0.
      tfaca(:) = 0.
      fwetl(:) = 0.
      fwets(:) = 0.
      fwetu(:) = 0.
      fwetlx(:) = 0.
      fwetsx(:) = 0.
      fwetux(:) = 0.
      rliql(:) = 0.
      rliqs(:) = 0.
      rliqu(:) = 0.

      allocate (var1(ndpts), var2(ndpts),&
               var3(ndpts), var4(ndpts),&
               var5(ndpts), var6(ndpts))
      var1(:) = 0.
      var2(:) = 0.
      var3(:) = 0.
      var4(:) = 0.
      var5(:) = 0.
      var6(:) = 0.

! Allocate variables in comcrop
      allocate 	(xirrig(lbeg:lend), xirriga(lbeg:lend), totirrig(lbeg:lend), df(lbeg:lend), vf(lbeg:lend), &
                cumvd(lbeg:lend), hdidx(lbeg:lend), gddfzcorn(lbeg:lend), gddfzsgc(lbeg:lend),             &
                gddfzsoy(lbeg:lend), gddfzwht(lbeg:lend), conspdate(lbeg:lend), conshybrid(lbeg:lend),     &
                cdays(lbeg:lend), ncyears(lbeg:lend), cropy(lbeg:lend), ik(lbeg:lend), gddpl15(lbeg:lend), &
                gddemerg(lbeg:lend), rm(lbeg:lend), consdays(lbeg:lend),              &
                maxcons(lbeg:lend), iniday(lbeg:lend), endday(lbeg:lend), gsdays(lbeg:lend))

      xirrig(:) = 0.
      xirriga(:) = 0.
      totirrig(:) = 0.
      df(:) = 0.
      vf(:) = 0.
      cumvd(:) = 0.
      hdidx(:) = 0.
      gddfzcorn(:) = 0.
      gddfzsgc(:) = 0.
      gddfzsoy(:) = 0.
      gddfzwht(:) = 0.
      conspdate(:) = 0.
      conshybrid(:) = 0.
      cdays(:) = 0.
      ncyears(:) = 0.
      cropy(:) = 0.
      ik(:) = 0.
      gddpl15(:) = 0.
      gddemerg(:) = 0.
      rm(:) = 0.
      consdays(:) = 0.
      maxcons(:) = 0.
      iniday(:) = 0
      endday(:) = 0
      gsdays(:) = 0

      allocate (aerial(npft), huileaf(npft), huigrain(npft))

      aerial(:) = 0.
      huileaf(:) = 0.
      huigrain(:) = 0.

      allocate (corndop(1,nrun+5), sgcdop(1,nrun+5), soydop(1,nrun+5), plmdop(1,nrun+5),   &
               whtdop(1,nrun+5), gddcorn(lbeg:lend,nrun+5), gddsgc(lbeg:lend,nrun+5),      &
               gddsgcp(lbeg:lend,2), gddsoy(lbeg:lend,nrun+5), gddwht(lbeg:lend,nrun+5),   &
               daygddc(1, 366), daygddsgc(1, 366), daygdds(1, 366))

      corndop(:,:) = 0.
      sgcdop(:,:) = 0.
      plmdop(:,:) = 0.
      soydop(:,:) = 0.
      whtdop(:,:) = 0.
      gddcorn(:,:) = 0.
      gddsgc(:,:) = 0.
      gddsgcp(:,:) = 0.
      gddsoy(:,:) = 0.
      gddwht(:,:) = 0.
      daygddc(:,:) = 0.
      daygddsgc(:,:) = 0.
      daygdds(:,:) = 0.


      allocate (cntops(lbeg:lend,npft),&
               cnrootvec(lbeg:lend,npft),&
               croplive(lbeg:lend,npft),&
               grnfraccrop(lbeg:lend,npft),&
	       gddplant(lbeg:lend,npft),&
	       gddtsoi(lbeg:lend,npft),&
	       gddmaturity(lbeg:lend,npft),&
	       thrlai(lbeg:lend,npft),&
	       peaklai(lbeg:lend,npft),&
	       hui(lbeg:lend,npft),&
	       phuf(lbeg:lend,npft),&
	       tlai(lbeg:lend,npft),&
	       templai(lbeg:lend,npft),&
	       harvidx(lbeg:lend,npft),&
	       fnleaf(lbeg:lend,npft),&
	       fnstem(lbeg:lend,npft),&
	       fnroot(lbeg:lend,npft),&
	       fngrain(lbeg:lend,npft),&
	       fnplant(lbeg:lend,npft),&
	       tnplant(lbeg:lend,npft),&
	       grainn(lbeg:lend,npft),&
	       cumlvs(lbeg:lend,npft),&
	       idpp(lbeg:lend,npft),&
	       dpgf(lbeg:lend,npft),&
	       cropyld(lbeg:lend,npft),&
	       dmleaf(lbeg:lend,npft),&
	       dmstem(lbeg:lend,npft),&
	       dmresidue(lbeg:lend,npft),&
	       dmyield(lbeg:lend,npft),&
	       dmcrop(lbeg:lend,npft),&
	       cropn(lbeg:lend, npft),&
	       cropfixn(lbeg:lend,npft),&
	       nconcl(lbeg:lend,npft),&
	       nconcs(lbeg:lend,npft),&
	       nconcr(lbeg:lend,npft),&
	       nconcg(lbeg:lend,npft),&
	       leafout(lbeg:lend,npft),&
	       cropplant(lbeg:lend,npft),&
	       croplaimx(lbeg:lend,npft),&
	       residuen(lbeg:lend, npft),&
	       dmroot(lbeg:lend,npft),&
	       hdate(lbeg:lend,npft),&
	       idppout(lbeg:lend,npft),&
	       pdate(lbeg:lend,npft),&
	       crmclim(lbeg:lend,npft),&
	       crmact(lbeg:lend,npft),&
	       crmplant(lbeg:lend,npft),&
	       ccdays(lbeg:lend,npft),&
	       grainday(lbeg:lend,npft),&
	       fertinput(lbeg:lend,npft),&
	       avehybrid(lbeg:lend,npft),&
	       pstart(lbeg:lend,npft),&
	       laidecl(lbeg:lend,npft),&
	       harvdate(lbeg:lend,npft),&
	       idop(lbeg:lend,npft),&
	       iavepdate(lbeg:lend,npft))

            cnrootvec(:,:) = cnroot
	    cntops(:,:) = 0.
	    croplive(:,:) = 0.
	    grnfraccrop(:,:) = 0.
	    gddplant(:,:) = 0.
	    gddtsoi(:,:) = 0.
	    gddmaturity(:,:) = 0.
	    thrlai(:,:) = 0.
	    peaklai(:,:) = 0.
	    hui(:,:) = 0.
	    phuf(:,:) = 0.
	    tlai(:,:) = 0.
	    templai(:,:) = 0.
	    harvidx(:,:) = 0.
	    fnleaf(:,:) = 0.
	    fnstem(:,:) = 0.
	    fnroot(:,:) = 0.
	    fngrain(:,:) = 0.
	    fnplant(:,:) = 0.
	    tnplant(:,:) = 0.
	    grainn(:,:) = 0.
!	    cumlvs(:,:) = 0.
	    idpp(:,:) = 0.
!	    dpgf(:,:) = 0.
	    cropyld(:,:) = 0.
	    dmleaf(:,:) = 0.
	    dmstem(:,:) = 0.
	    dmresidue(:,:) = 0.
	    dmyield(:,:) = 0.
	    dmcrop(:,:) = 0.
	    cropn(:,:) = 0.
	    cropfixn(:,:) = 0.
	    nconcl(:,:) = 0.
	    nconcs(:,:) = 0.
	    nconcr(:,:) = 0.
	    nconcg(:,:) = 0.
!	    leafout(:,:) = 0.
	    cropplant(:,:) = 0.
	    croplaimx(:,:) = 0.
	    residuen(:,:) = 0.
	    dmroot(:,:) = 0.
	    hdate(:,:) = 0.
	    idppout(:,:) = 0.
	    pdate(:,:) = 0.
	    crmclim(:,:) = 0.
	    crmact(:,:) = 0.
	    crmplant(:,:) = 0.
	    ccdays(:,:) = 0.
!	    grainday(:,:) = 0.
	    grainday(:,:) = 9999.
!	    fertinput(:,:) = 0.
	    avehybrid(:,:) = 0.
	    pstart(:,:) = 999.
	    laidecl(:,:) = 0.
	    harvdate(:,:) = 0.
	    idop(:,:) = 0.
	    iavepdate(:,:) = 0.

      allocate (cropout(lbeg:lend,npft,60))

               cropout(:,:,:) = 0.

      allocate (fleafwht(2),&
	       cfertmaize(2000),&
	       cfertsoy(2000),&
	       cfertsgc(2000),&
	       cfertwheat(2000))

	       fleafwht(:) = 0.
	       cfertmaize(:) = 0.
	       cfertsoy(:) = 0.
	       cfertsgc(:) = 0.
	       cfertwheat(:) = 0.

!Allocate variables in comnitr
      allocate (storedn(lbeg:lend), yrleach(lbeg:lend), aplantn(lbeg:lend),                                  &
               totimm(lbeg:lend), totmin(lbeg:lend), totnrel(lbeg:lend), ctot(lbeg:lend), ctoti(lbeg:lend),  &
               ftot(lbeg:lend), tsinp(lbeg:lend), tslay(lbeg:lend), taninp(lbeg:lend), tsnimm(lbeg:lend),    &
               tsnmob(lbeg:lend), dtnleach(lbeg:lend), dnileach(lbeg:lend), tpnuptake(lbeg:lend),            &
               totnvegn(lbeg:lend), drntot(lbeg:lend), ddrn(lbeg:lend), concn(lbeg:lend), assimn(lbeg:lend), &
               ydeposn(lbeg:lend), yfixsoin(lbeg:lend), yno3leach(lbeg:lend), snbalance(lbeg:lend),          &
               fixsoin(lbeg:lend))
      storedn(:) = 0.
      yrleach(:) = 0.
      aplantn(:) = 0.
      totimm(:) = 0.
      totmin(:) = 0.
      totnrel(:) = 0.
      ctot(:) = 0.
      ctoti(:) = 0.
      ftot(:) = 0.
      tsinp(:) = 0.
      tslay(:) = 0.
      taninp(:) = 0.
      tsnimm(:) = 0.
      tsnmob(:) = 0.
      dtnleach(:) = 0.
      dnileach(:) = 0.
      tpnuptake(:) = 0.
      totnvegn(:) = 0.
      drntot(:) = 0.
      ddrn(:) = 0.
      concn(:) = 0.
      assimn(:) = 0.
      ydeposn(:) = 0.
      yfixsoin(:) = 0.
      yno3leach(:) = 0.
      snbalance(:) = 0.
      fixsoin(:) = 0.

      allocate (totnuptake(lbeg:lend,npft), stressn(lbeg:lend,npft), totnfix(lbeg:lend,npft), &
               fixn(lbeg:lend,npft), fertnitro(lbeg:lend,npft))
      totnuptake(:,:) = 0.
      stressn(:,:) = 0.
      totnfix(:,:) = 0.
      fixn(:,:) = 0.
      fertnitro(:,:) = 0.

      allocate (daydrn(1,366), daynconc(1,366))

      daydrn = 0.
      daynconc = 0.

! Allocate variables in comveg
      allocate(q12(lbeg:lend),&
               q34(lbeg:lend),&
               sl(lbeg:lend),&
               ss(lbeg:lend),&
               su(lbeg:lend),&
               topparl(lbeg:lend),&
               topparu(lbeg:lend),&
               tl(lbeg:lend),&
               ts(lbeg:lend),&
               tu(lbeg:lend),&
               tlsub(lbeg:lend),&
               t12(lbeg:lend),&
               t34(lbeg:lend),&
               wliql(lbeg:lend),&
               wliqs(lbeg:lend),&
               wliqu(lbeg:lend),&
               wsnol(lbeg:lend),&
               wsnos(lbeg:lend),&
               wsnou(lbeg:lend),&
               disturbl(lbeg:lend))
               disturbl(:) = 0.

      allocate(agcub(lbeg:lend),&
               agcuc(lbeg:lend),&
               agcls(lbeg:lend),&
               agcl3(lbeg:lend),&
               agcl4(lbeg:lend),&
               ancub(lbeg:lend),&
               ancuc(lbeg:lend),&
               ancls(lbeg:lend),&
               ancl3(lbeg:lend),&
               ancl4(lbeg:lend),&
               totcondub(lbeg:lend),&
               totconduc(lbeg:lend),&
               totcondls(lbeg:lend),&
               totcondl3(lbeg:lend),&
               totcondl4(lbeg:lend))

      allocate(ciub(lbeg:lend),&
               ciuc(lbeg:lend),&
               cils(lbeg:lend),&
               cil3(lbeg:lend),&
               cil4(lbeg:lend),&
               csub(lbeg:lend),&
               csuc(lbeg:lend),&
               csls(lbeg:lend),&
               csl3(lbeg:lend),&
               csl4(lbeg:lend),&
               gsub(lbeg:lend),&
               gsuc(lbeg:lend),&
               gsls(lbeg:lend),&
               gsl3(lbeg:lend),&
               gsl4(lbeg:lend))

      allocate(agddl(lbeg:lend),&
               agddu(lbeg:lend),&
               fl(lbeg:lend),&
               fu(lbeg:lend),&
               gdd0(lbeg:lend),&
               daylength(lbeg:lend),&
               gdd0this(lbeg:lend),&
               gdd5(lbeg:lend),&
               gdd5this(lbeg:lend),&
               sapfrac(lbeg:lend),&
               tc(lbeg:lend),&
               tcthis(lbeg:lend),&
               tcmin(lbeg:lend),&
               totlail(lbeg:lend),&
               totlaiu(lbeg:lend),&
               totbiol(lbeg:lend),&
               totbiou(lbeg:lend),&
               tw(lbeg:lend),&
               twthis(lbeg:lend),&
               disturbf(lbeg:lend),&
               disturbo(lbeg:lend),&
               firefac(lbeg:lend),&
               tco2mic(lbeg:lend),&
               tco2root(lbeg:lend),&
               tneetot(lbeg:lend),&
               tnmin(lbeg:lend),&
               tnpptot(lbeg:lend),&
               tgpptot(lbeg:lend),&
               totalit(lbeg:lend),&
               totanlit(lbeg:lend),&
               totcmic(lbeg:lend),&
               totcsoi(lbeg:lend),&
               totfall(lbeg:lend),&
               totlit(lbeg:lend),&
               totnlit(lbeg:lend),&
               totnmic(lbeg:lend),&
               totnsoi(lbeg:lend),&
               totrlit(lbeg:lend),&
               totrnlit(lbeg:lend),&
               cleach(lbeg:lend),&
               tempu(lbeg:lend),&
               templ(lbeg:lend),&
               dropu(lbeg:lend),&
               dropls(lbeg:lend),&
               dropl4(lbeg:lend),&
               dropl3(lbeg:lend),&
               ancc4(lbeg:lend),&
               a3tdmin(lbeg:lend),&
               a11soiltd(lbeg:lend),&
               ancc3(lbeg:lend),&
               agcc3(lbeg:lend),&
               agcc4(lbeg:lend),&
               tsoiavg(lbeg:lend),&
               totcondc4(lbeg:lend),&
               totcondc3(lbeg:lend),&
               vegtype0(lbeg:lend))

      allocate(cbiol(lbeg:lend,npft),&
               cbior(lbeg:lend,npft),&
               cbiow(lbeg:lend,npft),&
               cbios(lbeg:lend,npft),&
               cbiog(lbeg:lend,npft),&
		       plaimx(lbeg:lend,npft))

      allocate(csoislo(lbeg:lend),&
               decompl(lbeg:lend),&
               decomps(lbeg:lend),&
               falll(lbeg:lend),&
               fallr(lbeg:lend),&
               fallw(lbeg:lend),&
               cdisturb(lbeg:lend),&
               greenfracl(lbeg:lend),&
               gdd12(lbeg:lend),&
               gdd11(lbeg:lend),&
               gdd10(lbeg:lend),&
               gdd8(lbeg:lend),&
               gdd0c(lbeg:lend),&
               gdd10this(lbeg:lend),&
               gdd8this(lbeg:lend),&
               gdd11this(lbeg:lend),&
               gdd12this(lbeg:lend),&
               gdd0cthis(lbeg:lend),&
               greenfracl3(lbeg:lend),&
               greenfracl4(lbeg:lend),&
               cic3(lbeg:lend),&
               cic4(lbeg:lend),&
               csc3(lbeg:lend),&
               csc4(lbeg:lend),&
               gsc3(lbeg:lend),&
               gsc4(lbeg:lend))

! FIXME: Useless. 0D model runs on only single point, No need to allocate it
!       as a vector but a single variable! Also, this seems to be used on a
!       single subroutine, so no need to declare on a module at all!
#ifdef SINGLE_POINT_MODEL
      allocate(cbioltot(npoi), cbiortot(npoi), cbiowtot(npoi))
#endif

      allocate(biomass(lbeg:lend,npft),&
               frac(lbeg:lend,npft),&
               plai(lbeg:lend,npft),&
               tnpp(lbeg:lend,npft),&
               tgpp(lbeg:lend,npft))

      allocate(lai(lbeg:lend,2),&
               sai(lbeg:lend,2),&
               zbot(lbeg:lend,2),&
               ztop(lbeg:lend,2),&
               ztopmx(lbeg:lend,2),&
               fallrsgc(lbeg:lend,3),&
		       falllsgc(lbeg:lend,3),&
		       htmx(lbeg:lend,2))

      ! dleaf and dstem are allocated on inland_prealloc
      allocate(orieh(2),oriev(2))

      allocate(froot(nsoilay,2))

      allocate(exist(lbeg:lend,npft))

      ! specla, aleaf, awood, aroot are allocated on inland_prealloc

      topparl(:) = 0.
      topparu(:) = 0.
      tl(:) = 0.
      t12(:) = 0.
      t34(:) = 0.

      agcub(:) = 0.
      agcuc(:) = 0.
      agcls(:) = 0.
      agcl3(:) = 0.
      agcl4(:) = 0.
      ancub(:) = 0.
      ancuc(:) = 0.
      ancls(:) = 0.
      ancl3(:) = 0.
      ancl4(:) = 0.
      ancc4(:) = 0.
      a3tdmin(:) = 0.
      a11soiltd(:) = 0.
      ancc3(:) = 0.
      agcc3(:) = 0.
      agcc4(:) = 0.
      tsoiavg(:) = 0.
      totcondc4(:) = 0.
      totcondc3(:) = 0.

      agddl(:) = 0.
      agddu(:) = 0.
      fl(:) = 0.
      fu(:) = 0.
      gdd0(:) = 0.
      daylength(:) = 0.
      gdd5(:) = 0.
      sapfrac(:) = 0.
      tc(:) = 0.
      tcmin(:) = 0.
      totlail(:) = 0.
      totlaiu(:) = 0.
      totbiol(:) = 0.
      totbiou(:) = 0.
      tw(:) = 0.
      disturbf(:) = 0.
      disturbo(:) = 0.
      firefac(:) = 0.
      tco2root(:) = 0.
      tgpptot(:) = 0.
      totnsoi(:) = 0.
      cleach(:) = 0.
      vegtype0(:) = 0.

      cbiol(:,:) = 0.
      cbior(:,:) = 0.
      cbiow(:,:) = 0.
      cbios(:,:) = 0.
      cbiog(:,:) = 0.

      csoislo(:) = 0.
      decompl(:) = 0.
      decomps(:) = 0.
      falll(:) = 0.
      fallr(:) = 0.
      fallw(:) = 0.

      biomass(:,:) = 0.
      frac(:,:) = 0.
      plai(:,:) = 0.
      tnpp(:,:) = 0.
      tgpp(:,:) = 0.

      lai(:,:) = 0.
      sai(:,:) = 0.
      zbot(:,:) = 0.
      ztop(:,:) = 0.
      ztopmx(:,:) = 0.
      fallrsgc(:,:) = 0.

      greenfracl(:) = 0.
      gdd12(:) = 0.
      gdd11(:) = 0.
      gdd10(:) = 0.
      gdd8(:) = 0.
      gdd0c(:) = 0.
      gdd10this(:) = 0.
      gdd8this(:) = 0.
      gdd11this(:) = 0.
      gdd12this(:) = 0.
      gdd0cthis(:) = 0.
      greenfracl3(:) = 0.
      greenfracl4(:) = 0.
      cic3(:) = 0.
      cic4(:) = 0.
      csc3(:) = 0.
      csc4(:) = 0.
      gsc3(:) = 0.
      gsc4(:) = 0.

      orieh(:) = 0.
      oriev(:) = 0.

      froot(:,:) = 0.

      exist(:,:) = 0.

! Allocate variables in comfire
      allocate(Pbio(lbeg:lend),&
               Pmoi(lbeg:lend),&
               Pign(lbeg:lend),&
               Pfire(lbeg:lend),&
               pfireyr(lbeg:lend),&
               srate(lbeg:lend),&
               lbratio(lbeg:lend),&
               abday(lbeg:lend),&
               burnfrac(lbeg:lend),&
               vegfrac(lbeg:lend))

      allocate(burnpft(lbeg:lend,npft))

      Pbio(:) = 0.
      Pmoi(:) = 0.
      Pign(:) = 0.
      Pfire(:) = 0.
      srate(:) = 0.
      lbratio(:) = 0.
      abday(:) = 0.
      burnfrac(:) = 0.
      vegfrac(:) = 0.
      pfireyr(:) = 0.
      burnpft(:,:) = 0.

! Allocate variables in comsoi
      allocate(wpud(lbeg:lend),&
               wipud(lbeg:lend),&
               z0soi(lbeg:lend),&
               albsav(lbeg:lend),&
               albsan(lbeg:lend),&
               stresstl(lbeg:lend),&
               stresstu(lbeg:lend),&
               heati(lbeg:lend),&
               heatg(lbeg:lend),&
               hvasug(lbeg:lend),&
               hvasui(lbeg:lend),&
               soihfl(lbeg:lend),&
               tg(lbeg:lend),&
               ti(lbeg:lend),&
               fwpudtot(lbeg:lend),&
               stre_tl(lbeg:lend),&
               stre_tu(lbeg:lend),&
               sice(lbeg:lend,nsoilay),&
               swater(lbeg:lend,nsoilay))

               wpud(:) = 0.
               wipud(:) = 0.
               z0soi(:) = 0.
               albsav(:) = 0.
               albsan(:) = 0.
               stresstl(:) = 0.
               stresstu(:) = 0.
               heati(:) = 0.
               heatg(:) = 0.
               hvasug(:) = 0.
               hvasui(:) = 0.
               soihfl(:) = 0.
               tg(:) = 0.
               ti(:) = 0.
               sice(:,:) = 0.
               swater(:,:) = 0.
! hsoi allocation moved to inland_prealloc.

      allocate(tsoi(lbeg:lend,nsoilay),&
               fracclay(lbeg:lend,nsoilay),&
               fracsand(lbeg:lend,nsoilay),&
               fracsilt(lbeg:lend,nsoilay),&
!               wsoi(lbeg:lend,nsoilay),&
               wisoi(lbeg:lend,nsoilay),&
               consoi(lbeg:lend,nsoilay),&
               csoi(lbeg:lend,nsoilay),&
               hydraul(lbeg:lend,nsoilay),&
               suction(lbeg:lend,nsoilay),&
               bex(lbeg:lend,nsoilay),&
               sfield(lbeg:lend,nsoilay),&
               swilt(lbeg:lend,nsoilay),&
               rhosoi(lbeg:lend,nsoilay),&
               poros(lbeg:lend,nsoilay),&
               porosflo(lbeg:lend,nsoilay),&
               stressl(lbeg:lend,nsoilay),&
               stressu(lbeg:lend,nsoilay),&
               upsoiu(lbeg:lend,nsoilay),&
               upsoil(lbeg:lend,nsoilay), &
	           cpwf(lbeg:lend,nsoilay),&
               tnuptake(lbeg:lend,nsoilay), &
               anuptake(lbeg:lend,nsoilay), &
               smsoil(lbeg:lend,nsoilay), &
               smsoln(lbeg:lend,nsoilay), &
               csoln(lbeg:lend,nsoilay), &
               fout(lbeg:lend,-nsoilay:nsoilay), &
               nout(lbeg:lend,-nsoilay:nsoilay))

      allocate(hflo(lbeg:lend,nsoilay+1),wflo(lbeg:lend,nsoilay+1))

      allocate(ibex(lbeg:lend,nsoilay))

      allocate(qglif(lbeg:lend,4))

      allocate(fwtop(lbeg:lend), fwpud(lbeg:lend),deposn(lbeg:lend))

      wflo(:,:) = 0.
      heati(:) = 0.
      heatg(:) = 0.
      hvasug(:) = 0.
      hvasui(:) = 0.
      soihfl(:) = 0.
      tg(:) = 0.
      fwtop(:) = 0.
      fwpud(:) = 0.
      deposn(:) = 0.

      tsoi(:,:) = 0.
      fracclay(:,:) = 0.
      fracsand(:,:) = 0.
      fracsilt(:,:) = 0.
!      wsoi(:,:) = 0.
      wisoi(:,:) = 0.
      consoi(:,:) = 0.
      porosflo(:,:) = 0.
      upsoiu(:,:) = 0.
      upsoil(:,:) = 0.
      tnuptake(:,:) = 0.
      anuptake(:,:) = 0.
      smsoil(:,:)=0.
      smsoln(:,:)=0.
      csoln(:,:)=0.
      fout(:,:)=0.
      nout(:,:)=0.
      hflo(:,:) = 0.

      qglif(:,:) = 0.

! Allocate variables in comsno
      allocate(fi(lbeg:lend),&
               tsno(lbeg:lend,nsnolay),&
               hsno(lbeg:lend,nsnolay))

      fi(:) = 0.
      tsno(:,:) = 0.
      hsno(:,:) = 0.

! Allocate variables in comatm
      allocate(coszen(lbeg:lend),&
               fira(lbeg:lend),&
               solsoi(lbeg:lend),&
               solad(lbeg:lend,nband),&
               solai(lbeg:lend,nband),&
               asurd(lbeg:lend,nband),&
               asuri(lbeg:lend,nband),&
               ua(lbeg:lend),&
               ux(lbeg:lend),&
               uy(lbeg:lend),&
               ta(lbeg:lend),&
               qa(lbeg:lend),&
               raina(lbeg:lend),&
               snowa(lbeg:lend),&
               psurf(lbeg:lend),&
               tmax(lbeg:lend),&
               tmin(lbeg:lend),&
               qd(lbeg:lend),&
               ud(lbeg:lend))

      coszen(:) = 0.
      fira(:) = 0.
!      solsoi(:) = 0.
      solad(:,:) = 0.
      solai(:,:) = 0.
      ua(:) = 0.
      ux(:) = 0.
      uy(:) = 0.
      ta(:) = 0.
      qa(:) = 0.
      raina(:) = 0.
      snowa(:) = 0.
      psurf(:) = 0.
      tmax(:) = 0.
      tmin(:) = 0.
      qd(:) = 0.
      if(irestart .eq. 1 .or. isimagro .eq.1)then
         ud(:) = 0.
      elseif(isimagro .eq. 0)then
         ud(:) = 0.
      endif

! comatm's weather generator specific variables
!---------------------------------------------------
!   There's no need to initialize these variables. This is done in the
! subroutine in a adequate frequency.
      allocate(precip(lbeg:lend),precipdaysum(lbeg:lend),iwet(lbeg:lend), &
               iwetdaysum(lbeg:lend),cloud(lbeg:lend))
      allocate(precipday(lbeg:lend,31),iwetday(lbeg:lend,31), &
               xstore(lbeg:lend,3))
      allocate(rh(lbeg:lend)) ! this variable is used only on diurnal subroutine
      xstore(:,:) = 0. !xstore needs to be initialized

! allocate seedvec, in comwork - values are set in main
      allocate(seedvec(lbeg:lend))
      seedvec(:) = 1

! Allocate variables in comhyd
      allocate(ginvap(lbeg:lend),gsuvap(lbeg:lend),gtrans(lbeg:lend),&
               gtransu(lbeg:lend),gtransl(lbeg:lend),grunof(lbeg:lend),&
               gdrain(lbeg:lend),gadjust(lbeg:lend))

      gtransu(:) = 0.
      gtransl(:) = 0.
      gadjust(:) = 0.

! combcs's weather generator specific variables
!-----------------------------------------------
! FIXME: MUST DECIDE::: in place of npoi, either lbeg:lend or mpt!!!
!   *** For CCSM3 convention, using lbeg:lend as the remaining variables at
!       combcs are using this label. This can yet impose problems to the model.
! Some 'missing' variables are used on readit thus allocated on inland_prealloc.
! Variables initialized on inland_prealloc will be associated to values on
! readit. Will these need to be zeroed?
      allocate(xinprecd(lbeg:lend),xintd(lbeg:lend),xintrngd(lbeg:lend), &
               xincldd(lbeg:lend),xinqd(lbeg:lend),xinwindd(lbeg:lend),&
               xintmaxd(lbeg:lend), xintmind(lbeg:lend))

! Allocate variables in comsum
      allocate(ndtime(numlv),nmtime(numlv),nytime(numlv))

      allocate(anytime(lbeg:lend))

      allocate(dwtot(lbeg:lend),wtotp(lbeg:lend),wtot(lbeg:lend))

      allocate(td(lbeg:lend),&
               adrain(lbeg:lend),&
               adsnow(lbeg:lend),&
               adaet(lbeg:lend),&
               adtrunoff(lbeg:lend),&
               adsrunoff(lbeg:lend),&
               addrainage(lbeg:lend),&
               adsnod(lbeg:lend),&
               adsnof(lbeg:lend),&
               adwsoi(lbeg:lend),&
               adwisoi(lbeg:lend),&
               adtsoi(lbeg:lend),&
               adwsoic(lbeg:lend),&
               adtsoic(lbeg:lend),&
               adco2mic(lbeg:lend),&
               adco2root(lbeg:lend),&
               adco2soi(lbeg:lend),&
               adco2ratio(lbeg:lend),&
               adnmintot(lbeg:lend),&
               adtlaysoi(lbeg:lend),&
               adwlaysoi(lbeg:lend),&
               adneetot(lbeg:lend),&
!               amtemp(lbeg:lend),&
!               amcloud(lbeg:lend),&
!               amrh(lbeg:lend),&
               amtrans(lbeg:lend),&
               amtratio(lbeg:lend),&
               amtotnleach(lbeg:lend),&
               amno3leach(lbeg:lend),&
!               amalbedo(lbeg:lend),&
               adtnpptot(lbeg:lend),&
               adpbio(lbeg:lend),&
               adpmoi(lbeg:lend),&
               adpign(lbeg:lend),&
               adpfire(lbeg:lend),&
               adsrate(lbeg:lend),&
               adabday(lbeg:lend),&
               adburnfrac(lbeg:lend),&
               adtrans(lbeg:lend),&
               adevap(lbeg:lend),&
               adtratio(lbeg:lend),&
!               adrh(lbeg:lend),&
!               adud(lbeg:lend),&
               adwsoi2(lbeg:lend))

      allocate(adburnpft(lbeg:lend,npft))

      allocate(amts2(lbeg:lend),&
               amrain(lbeg:lend),&
               amsnow(lbeg:lend),&
               amqa(lbeg:lend),&
               amaet(lbeg:lend),&
               amtransu(lbeg:lend),&
               amtransl(lbeg:lend),&
               amsuvap(lbeg:lend),&
               aminvap(lbeg:lend),&
               amtrunoff(lbeg:lend),&
               amsrunoff(lbeg:lend),&
               amdrainage(lbeg:lend),&
               amwsoi(lbeg:lend),&
               amwisoi(lbeg:lend),&
               amvwc(lbeg:lend),&
               amawc(lbeg:lend),&
               amtsoi(lbeg:lend),&
               amsnod(lbeg:lend),&
               amsnof(lbeg:lend),&
               amlaiu(lbeg:lend),&
               amlail(lbeg:lend),&
               amsolar(lbeg:lend),&
               amreflect(lbeg:lend),&
               amirdown(lbeg:lend),&
               amirup(lbeg:lend),&
               amsens(lbeg:lend),&
               amlatent(lbeg:lend),&
               amnpptot(lbeg:lend),&
               amneetot(lbeg:lend),&
               amco2mic(lbeg:lend),&
               amco2root(lbeg:lend),&
               amco2soi(lbeg:lend),&
               amco2ratio(lbeg:lend),&
               amnmintot(lbeg:lend),&
               ampbio(lbeg:lend),&
               ampmoi(lbeg:lend),&
               ampign(lbeg:lend),&
               ampfire(lbeg:lend),&
               amsrate(lbeg:lend),&
               amabday(lbeg:lend),&
               amburnfrac(lbeg:lend),&
               a10tmin(lbeg:lend),&
               a5tmin(lbeg:lend),&
               a5td(lbeg:lend))

      allocate(adcsoln(lbeg:lend,nsoilay),&
               adwisoilay(lbeg:lend,nsoilay),&
               adwsoilay(lbeg:lend,nsoilay))
!               adtsoilay(lbeg:lend,nsoilay),&
!               adupsoil(lbeg:lend,nsoilay))

      allocate(amnpp(lbeg:lend,npft), &
               amburnpft(lbeg:lend,npft),&
               aybprod(lbeg:lend,npft), &
               ayabprod(lbeg:lend,npft), &
	       ayrprod(lbeg:lend,npft), &
	       aylprod(lbeg:lend,npft))

      aybprod(:,:) = 0.
      ayrprod(:,:) = 0.
      aylprod(:,:) = 0.

      allocate(ayprcp(lbeg:lend),&
               ayaet(lbeg:lend),&
               aytrans(lbeg:lend),&
               aytrunoff(lbeg:lend),&
               aysrunoff(lbeg:lend),&
               aydrainage(lbeg:lend),&
               aywsoi(lbeg:lend),&
               aywisoi(lbeg:lend),&
               ayvwc(lbeg:lend),&
               ayawc(lbeg:lend),&
               aytsoi(lbeg:lend),&
               ayrratio(lbeg:lend),&
               aytratio(lbeg:lend),&
               aysolar(lbeg:lend),&
               ayreflect(lbeg:lend),&
               ayirdown(lbeg:lend),&
               ayirup(lbeg:lend),&
               aysens(lbeg:lend),&
               aylatent(lbeg:lend),&
               aystresstu(lbeg:lend),&
               aystresstl(lbeg:lend),&
               ayanpptot(lbeg:lend),&
               aynpptot(lbeg:lend),&
               aygpptot(lbeg:lend),&
               ayalit(lbeg:lend),&
               ayblit(lbeg:lend),&
               aycsoi(lbeg:lend),&
               aycmic(lbeg:lend),&
               ayanlit(lbeg:lend),&
               aybnlit(lbeg:lend),&
               aynsoi(lbeg:lend),&
               ayneetot(lbeg:lend),&
               ayco2mic(lbeg:lend),&
               ayco2root(lbeg:lend),&
               ayco2soi(lbeg:lend),&
               aynmintot(lbeg:lend),&
               ayrootbio(lbeg:lend),&
               aymintot(lbeg:lend),&
               ayimmtot(lbeg:lend),&
               aynreltot(lbeg:lend),&
               ayalbedo(lbeg:lend),&
               caccount(lbeg:lend),&
               aypbio(lbeg:lend),&
               aypmoi(lbeg:lend),&
               aypign(lbeg:lend),&
               aypfire(lbeg:lend),&
               aysrate(lbeg:lend),&
               ayabday(lbeg:lend),&
               ayburnfrac(lbeg:lend))

      allocate(ayanpp(lbeg:lend,npft),&
               aynpp(lbeg:lend,npft),&
               aygpp(lbeg:lend,npft),&
               ayburnpft(lbeg:lend,npft),&
               adnpp(lbeg:lend,npft))

      allocate(a10td(lbeg:lend),&
               a10ancub(lbeg:lend),&
               a10ancuc(lbeg:lend),&
               a10ancls(lbeg:lend),&
               a10ancl3(lbeg:lend),&
               a10ancl4(lbeg:lend),&
               a10scalparamu(lbeg:lend),&
               a10scalparaml(lbeg:lend),&
               a10daylightu(lbeg:lend),&
               a10ts(lbeg:lend),&
               a10ancc3(lbeg:lend),&
               a10ancc4(lbeg:lend),&
               a10daylightl(lbeg:lend))

      allocate(ynleach(lbeg:lend))

      ndtime(:) = 0
      nmtime(:) = 0
      nytime(:) = 0

      anytime(:) = 0.

! moved to inland_parameters (subroutine) allocation area on this file
!     garea(:) = 0.

      dwtot(:) = 0.
      wtotp(:) = 0.
      wtot(:) = 0.

      adtsoic(:) = 0.
      adco2ratio(:) = 0.
      adneetot(:) = 0.
      adtnpptot(:) = 0.
      adpbio(:) = 0.
      adpmoi(:) = 0.
      adpign(:) = 0.
      adpfire(:) = 0.
      adsrate(:) = 0.
      adabday(:) = 0.
      adburnfrac(:) = 0.
      adburnpft(:,:) = 0.
!      amtemp(:) = 0.
!      amcloud(:) = 0.
!      amrh(:) = 0.
      amtrans(:) = 0.
      amtratio(:) = 0.
      amtotnleach(:) = 0.
      amno3leach(:) = 0.
!      amalbedo(:) = 0.
      adtnpptot(:) = 0.
      adtrans(:) = 0.
      adevap(:) = 0.
      adtratio(:) = 0.
!      adrh(:) = 0.
!      adud(:) = 0.
      adwsoi2(:) = 0.
      amts2(:) = 0.
      amnpptot(:) = 0.
      amneetot(:) = 0.
      amco2soi(:) = 0.
      amco2ratio(:) = 0.
      ampbio(:) = 0.
      ampmoi(:) = 0.
      ampign(:) = 0.
      ampfire(:) = 0.
      amsrate(:) = 0.
      amabday(:) = 0.
      amburnfrac(:) = 0.
      amburnpft(:,:) = 0.
      a10tmin(:) = 0.
      a5tmin(:) = 0.
      a5td(:) = 0.

      adcsoln(:,:) = 0.
      adwisoilay(:,:) = 0.
      adwsoilay(:,:) = 0.
!      adtsoilay(:,:) = 0.
!      adupsoil(:,:) = 0.

      ayprcp(:) = 0.
      ayaet(:) = 0.
      aytrans(:) = 0.
      aytrunoff(:) = 0.
      aysrunoff(:) = 0.
      aydrainage(:) = 0.
      aywsoi(:) = 0.
      aywisoi(:) = 0.
      ayvwc(:) = 0.
      ayawc(:) = 0.
      aytsoi(:) = 0.
      ayrratio(:) = 0.
      aytratio(:) = 0.
      aysolar(:) = 0.
      ayreflect(:) = 0.
      ayirdown(:) = 0.
      ayirup(:) = 0.
      aysens(:) = 0.
      aylatent(:) = 0.
      aystresstu(:) = 0.
      aystresstl(:) = 0.
      ayanpptot(:) = 0.
      aynpptot(:) = 0.
      aygpptot(:) = 0.
      ayalit(:) = 0.
      ayblit(:) = 0.
      aycsoi(:) = 0.
      aycmic(:) = 0.
      ayanlit(:) = 0.
      aybnlit(:) = 0.
      aynsoi(:) = 0.
      ayneetot(:) = 0.
      ayco2mic(:) = 0.
      ayco2root(:) = 0.
      ayco2soi(:) = 0.
      aynmintot(:) = 0.
      ayrootbio(:) = 0.
      aymintot(:) = 0.
      ayimmtot(:) = 0.
      aynreltot(:) = 0.
      ayalbedo(:) = 0.
      caccount(:) = 0.
      aypbio(:) = 0.
      aypmoi(:) = 0.
      aypign(:) = 0.
      aypfire(:) = 0.
      aysrate(:) = 0.
      ayabday(:) = 0.
      ayburnfrac(:) = 0.
      ayburnpft(:,:) = 0.

      ayanpp(:,:) = 0.
      aynpp(:,:) = 0.
      aygpp(:,:) = 0.
      adnpp(:,:) = 0.

      a10td(:) = 0.
      a10ancub(:) = 0.
      a10ancuc(:) = 0.
      a10ancls(:) = 0.
      a10ancl3(:) = 0.
      a10ancl4(:) = 0.
      a10scalparamu(:) = 0.
      a10scalparaml(:) = 0.
      a10daylightu(:) = 0.
      a10ts(:) = 0.
      a10ancc3(:) = 0.
      a10ancc4(:) = 0.
      a10daylightl(:) = 0.

      ynleach(:) = 0.

! All comtex module variables is allocated on inland_prealloc!

! Allocate variables in compft
!   vmax_pft, tauleaf, tauroot, tauwood0, TminL, TminU, Twarm, GDD, plai_init
! are all allocated on inland_prealloc
      allocate(tauwood(npoi, npft))
      tauwood(:,:) = 0.

#ifdef SINGLE_POINT_MODEL

! Allocate variables in comforc (forcings for the single point model)
      allocate(xta(npoi,dimforc),xprec(npoi,dimforc),xua(npoi,dimforc),  &
               xcld(npoi,dimforc),xqa(npoi,dimforc),                     &
               xsmoi(npoi,nsoilay,dimforc),xstemp(npoi,nsoilay,dimforc), &
               xsin(npoi,dimforc),xlin(npoi,dimforc),xlati(npoi))

! TODO: zero only variables that need to be zeroed. - fzm
      xta(:,:) = 0.
      xprec(:,:) = 0.
      xua(:,:) = 0.
      xcld(:,:) = 0.
      xqa(:,:) = 0.
      xsmoi(:,:,:) = 0.
      xstemp(:,:,:) = 0.
      xsin(:,:) = 0.
      xlin(:,:) = 0.
      xlati(:) = 0.
#endif /* SINGLE_POINT_MODEL */

      return
end subroutine alloc
