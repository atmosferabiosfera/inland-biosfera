#include "inland_config.h"
subroutine wdiag(i, istep)
! ---------------------------------------------------------------------
! CALLED by main
! CALLS nothing
!
! This subroutine writes diagnostic output to the "diag" files
      use inland_parameters
      use inland_control, only: iyear, imonth, iday
      use inland_com1d
      use inland_comveg
      use inland_comhyd
      use inland_comatm
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comdiag

      real*8      diag(nvars)              ! diagnostic variables
      real*8      rowdata(16)              ! rows of data for output

      integer   i,       & ! # diagnostic file requested
                ivar,    & ! index
                kolumn,  & ! columns for output
                countvar   !

      diag(1) = adco2mic(ndiagpt(i))
      diag(2) = adco2root(ndiagpt(i))
      diag(3) = adco2soi(ndiagpt(i))
      diag(4) = adneetot(ndiagpt(i))
      diag(5) = adtsoic(ndiagpt(i))
      diag(6) = adwsoic(ndiagpt(i))
      diag(7) = amco2root(ndiagpt(i))
      diag(8) = amco2soi(ndiagpt(i))
      diag(9) = ancl3(ndiagpt(i))
      diag(10) = ancl4(ndiagpt(i))
      diag(11) = ancub(ndiagpt(i))
      diag(12) = ancuc(ndiagpt(i))
      diag(13) = asurd(ndiagpt(i),1)
      diag(14) = asurd(ndiagpt(i),2)
      diag(15) = asuri(ndiagpt(i),1)
      diag(16) = asuri(ndiagpt(i),2)
      diag(17) = ayco2soi(ndiagpt(i))
      diag(18) = biomass(ndiagpt(i),1)
      diag(19) = biomass(ndiagpt(i),2)
      diag(20) = biomass(ndiagpt(i),3)
      diag(21) = biomass(ndiagpt(i),4)
      diag(22) = biomass(ndiagpt(i),5)
      diag(23) = biomass(ndiagpt(i),6)
      diag(24) = biomass(ndiagpt(i),7)
      diag(25) = biomass(ndiagpt(i),8)
      diag(26) = biomass(ndiagpt(i),9)
      diag(27) = biomass(ndiagpt(i),10)
      diag(28) = biomass(ndiagpt(i),11)
      diag(29) = biomass(ndiagpt(i),12)
      diag(30) = cloud(ndiagpt(i))
      diag(31) = coszen(ndiagpt(i))
      diag(32) = fi(ndiagpt(i))
      diag(33) = fira(ndiagpt(i))
      diag(34) = firb(ndiagpt(i))
      diag(35) = frac(ndiagpt(i),1)
      diag(36) = frac(ndiagpt(i),2)
      diag(37) = frac(ndiagpt(i),3)
      diag(38) = frac(ndiagpt(i),4)
      diag(39) = frac(ndiagpt(i),5)
      diag(40) = frac(ndiagpt(i),6)
      diag(41) = frac(ndiagpt(i),7)
      diag(42) = frac(ndiagpt(i),8)
      diag(43) = frac(ndiagpt(i),9)
      diag(44) = frac(ndiagpt(i),10)
      diag(45) = frac(ndiagpt(i),11)
      diag(46) = frac(ndiagpt(i),12)
      diag(47) = gadjust(ndiagpt(i))
      diag(48) = gdrain(ndiagpt(i))
      diag(49) = ginvap(ndiagpt(i))
      diag(50) = grunof(ndiagpt(i))
      diag(51) = gsuvap(ndiagpt(i))
      diag(52) = gtrans(ndiagpt(i))
      diag(53) = gtransl(ndiagpt(i))
      diag(54) = gtransu(ndiagpt(i))
      diag(55) = hsno(ndiagpt(i),1)
      diag(56) = hsno(ndiagpt(i),2)
      diag(57) = hsno(ndiagpt(i),3)
      diag(58) = hsnotop
      diag(59) = plai(ndiagpt(i),1)
      diag(60) = plai(ndiagpt(i),2)
      diag(61) = plai(ndiagpt(i),3)
      diag(62) = plai(ndiagpt(i),4)
      diag(63) = plai(ndiagpt(i),5)
      diag(64) = plai(ndiagpt(i),6)
      diag(65) = plai(ndiagpt(i),7)
      diag(66) = plai(ndiagpt(i),8)
      diag(67) = plai(ndiagpt(i),9)
      diag(68) = plai(ndiagpt(i),10)
      diag(69) = plai(ndiagpt(i),11)
      diag(70) = plai(ndiagpt(i),12)
      diag(71) = precip(ndiagpt(i))
      diag(72) = psurf(ndiagpt(i))
      diag(73) = qa(ndiagpt(i))
      diag(74) = raina(ndiagpt(i))
      diag(75) = snowa(ndiagpt(i))
      diag(76) = solad(ndiagpt(i),1)
      diag(77) = solad(ndiagpt(i),2)
      diag(78) = solai(ndiagpt(i),1)
      diag(79) = solai(ndiagpt(i),2)
      diag(80) = ta(ndiagpt(i))
      diag(81) = td(ndiagpt(i))
      diag(82) = tl(ndiagpt(i))
      diag(83) = totalit(ndiagpt(i))
      diag(84) = totcmic(ndiagpt(i))
      diag(85) = totcsoi(ndiagpt(i))
      diag(86) = totrlit(ndiagpt(i))
      diag(87) = tmax(ndiagpt(i))
      diag(88) = ts(ndiagpt(i))
      diag(89) = tsno(ndiagpt(i),1)
      diag(90) = tsno(ndiagpt(i),2)
      diag(91) = tsno(ndiagpt(i),3)
      diag(92) = tsoi(ndiagpt(i),1)
      diag(93) = tsoi(ndiagpt(i),2)
      diag(94) = tsoi(ndiagpt(i),3)
      diag(95) = tsoi(ndiagpt(i),4)
      diag(96) = tsoi(ndiagpt(i),5)
      diag(97) = tsoi(ndiagpt(i),6)
      diag(98) = tu(ndiagpt(i))
      diag(99) = iwet(ndiagpt(i))
      diag(100) = wipud(ndiagpt(i))
      diag(101) = wisoi(ndiagpt(i),1)
      diag(102) = wisoi(ndiagpt(i),2)
      diag(103) = wisoi(ndiagpt(i),3)
      diag(104) = wisoi(ndiagpt(i),4)
      diag(105) = wisoi(ndiagpt(i),5)
      diag(106) = wisoi(ndiagpt(i),6)
      diag(107) = wliql(ndiagpt(i))
      diag(108) = wliqs(ndiagpt(i))
      diag(109) = wliqu(ndiagpt(i))
      diag(110) = wpud(ndiagpt(i))
      diag(111) = wsnol(ndiagpt(i))
      diag(112) = wsnos(ndiagpt(i))
      diag(113) = wsnou(ndiagpt(i))
      diag(114) = wsoi(ndiagpt(i),1)
      diag(115) = wsoi(ndiagpt(i),2)
      diag(116) = wsoi(ndiagpt(i),3)
      diag(117) = wsoi(ndiagpt(i),4)
      diag(118) = wsoi(ndiagpt(i),5)
      diag(119) = wsoi(ndiagpt(i),6)
      diag(120) = biomass(ndiagpt(i),13)
      diag(121) = biomass(ndiagpt(i),14)
      diag(122) = biomass(ndiagpt(i),15)
      diag(123) = frac(ndiagpt(i),13)
      diag(124) = frac(ndiagpt(i),14)
      diag(125) = frac(ndiagpt(i),15)
      diag(126) = plai(ndiagpt(i),13)
      diag(127) = plai(ndiagpt(i),14)
      diag(128) = plai(ndiagpt(i),15)

      rowdata(1) = float(iyear)
      rowdata(2) = float(imonth)
      rowdata(3) = float(iday)
      rowdata(4) = float(istep)

      kolumn = 4
      countvar = 0
      do 30 ivar = 1, nvars
         if( ldiag(ivar,i) .eq. 1) then
            countvar = countvar + 1
            kolumn = kolumn + 1
            rowdata(kolumn) = diag(ivar)
         endif
30    continue

      if (countvar .lt. 12) then
         do 40 n = countvar + 5, 16
            rowdata(n) = -999.
40       continue
      endif
      write(i+20, 3000) rowdata
3000  format(1x,f5.0,2x,f3.0,3x,f7.0,1x,f5.0,1x,12(e12.3,1x))
      return
end subroutine wdiag
