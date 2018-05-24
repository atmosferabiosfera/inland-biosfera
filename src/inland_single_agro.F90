#include "inland_config.h"

subroutine single_agro

      use inland_comatm, only: fira, solai, asuri, solad, asurd, ta, solsoi
      use inland_com1d, only: firb, fsena, fvapa, abupi, abupd, flx
      use inland_parameters, only: dtime, hvap
      use inland_comveg, only: use, topparu, topparl, tneetot, plai, tl
      use inland_control, only: istep, iyear, jday
      use inland_comsoi, only: soihfl, tsoi
      use inland_comhyd, only: gsuvap, gtrans
      use inland_comcrop, only: grnfraccrop
      use inland_comhour, only:imetyear

      implicit none

      real *8 rn,          & ! net radiation flux (short and long waves)
              pari,        & ! incoming PAR
              apar,        & ! APAR
              paro,        & ! outgoing PAR
              swin,        & ! incoming solar radiation (W/m²)
              swout          ! reflect solar radiation (W/m²)

           rn = solad(1,1) * (1 - asurd(1,1)) + solad(1,2) * (1 - asurd(1,2)) + &
                solai(1,1) * (1 - asuri(1,1)) + solai(1,2) * (1 - asuri(1,2)) + & 
                fira(1) - firb(1)

           swin = solad(1,1) + solad(1,2) + solai(1,1) + solai(1,2)
    
           swout= solad(1,1) * asurd(1,1) + solad(1,2) * asurd(1,2)  + solai(1,1) * &
                  asuri(1,1) + solai(1,2) * asuri(1,2)

           pari = (solad(1,1) + solai(1,1)) * 4.59
           apar = (topparl(1) + topparu(1)) * 4.59
	              
           paro = (solad(1,1) * asurd(1,1) + solai(1,1) * asuri(1,1))* 4.59 

!           if((iyear.eq.2009.and.jday.gt.300).or.(iyear.eq.2010.and.jday.le.120)) then
           if(imetyear .ne. 9999) then

                write(43,9200) iyear, jday, istep-1, swin, swout, pari, paro,  &
                               apar, rn,-fsena(1),-fvapa(1)*hvap,-soihfl(1),   & 
                               -tneetot(1)*1e6,plai(1,16)  
         
                if(istep-1.ge.8.and.istep-1.le.18) then  !diurnal
             
                     write(229,*)iyear,jday,istep-1,-fsena(1),- fvapa(1)*hvap, &
                                 -tneetot(1)*1e6,ta(1)-273.16,tl(1)-273.16

	             use(4,istep-1)=-tneetot(1)*1e6
             
                elseif(istep-1.eq.19) then
	
                     write(27,173)jday,iyear,plai(1,13)*grnfraccrop(1,13),                    &
                                  use(4,8),use(4,9),use(4,10),use(4,11),use(4,12),            &
                                  use(4,13),use(4,14),use(4,15),use(4,16),use(4,17),use(4,18) 

                     173   format (1x,i3,1x,i4,1x,f4.2,11(1x,f6.2))	

                endif
           endif

9200   format (3(i4,'  '),10(f9.3,'  '),f3.1)  


             flx(1)= flx(1) + swin   *(dtime/10**6) 
             flx(2)= flx(2) + swout  *(dtime/10**6)
             flx(3)= flx(3) + pari   
             flx(4)= flx(4) + paro   
             flx(5)= flx(5) + apar   
             flx(6)= flx(6) + rn   *(dtime/10**6)

             flx(7)= flx(7) - fsena(1)  *(dtime/10**6)     

             flx(8)= flx(8) - fvapa(1)*hvap *(dtime/10**6)

             flx(9) = flx(9) - soihfl(1)*(dtime/10**6)
		
             flx(10) = flx(10) -tneetot(1)*1e6 

             flx(11) = flx(11)+ solsoi(1)   *(dtime/10**6)

             flx(12) = flx(12)+ gsuvap(1) * dtime
             flx(13) = flx(13)+ gtrans(1) * dtime

             flx(14) = flx(14)+ tsoi(1,1) !-273.18

             flx(15) = flx(15)+ ((tsoi(1,4)+tsoi(1,5))/2.)!-273.18 

      return
end subroutine single_agro
