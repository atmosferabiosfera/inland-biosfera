#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine inisum (irestart)
! ---------------------------------------------------------------------
! does initialization for time averaging
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comveg

      implicit none
! ---------------------------------------------------------------------
! input argument
      integer irestart     ! 0: not restart of INLAND, 1 restart of INLAND

! local variables
      integer i, k         ! loop indices

      if (irestart .eq. 0) then

! initialize total water content in soil+snow+vegetation
! (for mass conservation check)
         do 20 i = lbeg, lend
            wtot(i) = (wliqu(i)+wsnou(i)) * fu(i) * 2.0 * lai(i,2) + &
                      (wliqs(i)+wsnos(i)) * fu(i) * 2.0 * sai(i,2) + &
                      (wliql(i)+wsnol(i)) * fl(i) * 2.0 * &
                      (lai(i,1) + sai(i,1)) * (1. - fi(i))
            wtot(i) = wtot(i) + wpud(i) + wipud(i)
            do 10 k = 1, nsoilay
               wtot(i) = wtot(i) + &
                         poros(i,k)*wsoi(i,k)*(1.-wisoi(i,k))*hsoi(k)*rhow + &
                         poros(i,k)*wisoi(i,k)*hsoi(k)*rhow
10          continue
            do 21 k = 1, nsnolay
               wtot(i) = wtot(i) + fi(i)*rhos*hsno(i,k)
21          continue
20       continue
      endif

! Daily means
      adrain(:)=0.0
      adsnow(:)=0.0
      adaet(:)=0.0
      adtrunoff(:)=0.0
      adsrunoff(:)=0.0
      addrainage(:)=0.0
      adsnod(:)=0.0
      adsnof(:)=0.0
      adwsoi(:)=0.0
      adtsoi(:)=0.0
      adwisoi(:)=0.0
      adtlaysoi(:)=0.0
      adwlaysoi(:)=0.0
      adwsoic(:)=0.0
      adtsoic(:)=0.0
      adco2mic(:)=0.0
      adco2root(:)=0.0
      adco2soi(:)=0.0
      adnmintot(:)=0.0
      if (irestart .eq. 0) then
         decompl(:)=0.0
         decomps(:)=0.0
      end if

! Monthly mean quanties
      amrain(:)=0.0
      amsnow(:)=0.0
      amaet(:)=0.0
      amtrunoff(:)=0.0
      amsrunoff(:)=0.0
      amdrainage(:)=0.0
      amqa(:)=0.0
      amsolar(:)=0.0
      amirup(:)=0.0
      amirdown(:)=0.0
      amreflect(:)=0.0
      amsens(:)=0.0
      amlatent(:)=0.0
      amtransu(:)=0.0
      amtransl(:)=0.0
      amsuvap(:)=0.0
      aminvap(:)=0.0
      amlaiu(:)=0.0
      amlail(:)=0.0
      amtsoi(:)=0.0
      amwsoi(:)=0.0
      amwisoi(:)=0.0
      amvwc(:)=0.0
      amawc(:)=0.0
      amsnod(:)=0.0
      amsnof(:)=0.0
      amco2mic(:)=0.0
      amco2root(:)=0.0
      amnmintot(:)=0.0
      amnpp(:,:)=0.0

! Annual mean quantities
      if (irestart .eq. 0) then
         aysolar(:)=0.0
         ayreflect(:)=0.0
         ayirup(:)=0.0
         ayirdown(:)=0.0
         aysens(:)=0.0
         aylatent(:)=0.0
         ayprcp(:)=0.0
         ayaet(:)=0.0
         aytrans(:)=0.0
         aysrunoff(:)=0.0
         aydrainage(:)=0.0
         dwtot(:)=0.0
         ayco2mic(:)=0.0
         ayco2root(:)=0.0
         ayalit(:)=0.0
         ayblit(:)=0.0
         aycsoi(:)=0.0
         ayanlit(:)=0.0
         aybnlit(:)=0.0
         aynsoi(:)=0.0
         aygpp(:,:)=0.0
         aynpp(:,:)=0.0
         aywsoi(:)=0.0
         aywisoi(:)=0.0
         aytsoi(:)=0.0
         ayvwc(:)=0.0
         ayawc(:)=0.0
!         aystresstu(:)=0.0
!         aystresstl(:)=0.0
         firefac(:)=0.0
         ayrootbio(:)=0.0
         aynmintot(:)=0.0
         aycmic(:)=0.0
      endif
      return
end subroutine inisum 
