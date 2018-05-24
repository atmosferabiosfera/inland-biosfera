#include "inland_config.h"

#ifndef SINGLE_POINT_MODEL
#error "This subroutine should ONLY be compiled for 0D INLAND model option."
#endif

subroutine readforc (isoilforc)
      use inland_com1d, only: za
      use inland_parameters, only: dtime, ndaypm, nsoilay, pi
      use inland_control, only: datadir, indir
      use inland_control, only: imonth0, iday0 
      use inland_combcs, only: xintopo, xinveg, deltat
      use inland_comsoi, only: sand, clay, tsoi
      use inland_comveg, only: exist, tc, tw, gdd0, gdd5, tcmin, ztop, zbot
      use inland_comforc, only: xlati, xua, xta, xprec, xsin, xlin, xqa, dimforc

      implicit none

! Arguments
      integer  isoilforc!, & ! soil physics dynamics (1) or static (0)
!               dimforc      ! dimension (number of lines) of input data

! local variables
      integer  it0,       & ! define the first step of a full data month
               it, month, ij, niter

!  MHC on Dec-6-2004: xtd should be defined at least as long as the number 
!  of days
      real*8   xtd(365),  & ! daily temperature over the first 365 days
               zsumtd, xtm
      integer dum1, dum2, dum3

      character*21 dum

! Reading soil and vegetation characteristics, and initial conditions
      niter = int (86400./dtime)

      open (10,status='unknown',file=trim(indir)//'/params/single_point_parameters')
      read(10,*) imonth0
      read(10,*) iday0
      read(10,*) xintopo(1)     !necessary to calculate sfc pressure
      read(10,*) xlati(1)
      read(10,*) xinveg(1)
      read(10,*) za(1)
      read(10,*) tsoi(1,1)
      read(10,*) ztop(1,2)
      read(10,*) zbot(1,2)
      read(10,*) exist(1,:)

! convert xlat in radians
!   With the calculation here (not in diurnal), the xlati coefficient is not re-
! converted systematically, thus turning into 0.
      xlati(1) = xlati(1) * pi / 180.0
 
!      if (isimveg.eq.0) then
!         read(10,*) (laium(1,k),k=1,12)
!         read(10,*) (fum(1,k),k=1,12)
!      else
!         read(10,*)             !skip 2 lines
!         read(10,*)
!      endif
      read(10,*) (sand(1,ij),ij=1,nsoilay)
      read(10,*) (clay(1,ij),ij=1,nsoilay)

      close(10)

      do 899 it=2,nsoilay

        tsoi(1,it)=tsoi(1,1)

899 continue      


      write(*,*) 'soil and vegetation general info read'

! Reading atmospheric forcing
      open (20,status='old',file=trim(datadir)//'/clim-input')

      do it=1, dimforc

         if (isoilforc.eq.1) then
            read(20,10001,end=900) dum1,dum2,dum3,xua(1,it),xta(1,it),       &
                                   xprec(1,it),xsin(1,it),xlin(1,it), &
                                   xqa(1,it)
         else
            read(20,10001,end=900) dum1,dum2,dum3,xua(1,it),xta(1,it),       &
                                   xprec(1,it),xsin(1,it),xlin(1,it), &
                                   xqa(1,it)
         endif
      enddo
10001 format(3i8,6f8.3)

900    close(20)
      write(*,*) 'NOTICE: ',dimforc/24,' days of climate and soil data read from clim-input.'

! Calculation of some variable to shortcut weather.f
! Calculation of td over the first 365 days of simulation
      zsumtd = 0
      ij = 0
      do 901 it = 1, dimforc
         zsumtd = zsumtd + xta(1,it)/float(niter)
         if (mod(it,niter) .eq. 0) then
            if (ij.ge.size(xtd)) goto 901 ! do nothing if already have the desired days.
            ij = ij + 1
            xtd(ij) = zsumtd
            zsumtd = 0.
         endif
901   continue

! calculation of gdd5, gdd0, tc, tw, and tcmin
      it0 = ndaypm(imonth0) - iday0 + 1
      tc(1) = 100.
      tw(1) = -50.
      gdd0(1) = 0.
      gdd5(1) = 0.
      do 902 it = 1, 12
         month = imonth0+it
         if (month.gt.12) month = month - 12
         do 903 ij = it0, it0 + ndaypm(month)
            if (ij.gt.size(xtd)) goto 903 ! do nothing beyond xtd size
            xtm = xtm + xtd(ij) / ndaypm(month)
            gdd0(1) = gdd0(1) + max(dble(0.0),xtd(ij))
            gdd5(1) = gdd5(1) + max(dble(0.0),xtd(ij))
903      continue
         it0=it0 + ndaypm(month)
         tc(1) = min(tc(1), xtm)
         tw(1) = max(tw(1), xtm)
902   continue

      return
end subroutine readforc
