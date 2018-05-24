#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine summonth (mcsec, loopi, kpti, kptj)
! ---------------------------------------------------------------------
!
! first convert to units that make sense for output
!
!   - convert all temperatures to deg c
!   - convert all liquid or vapor fluxes to mm/day
!   - redefine upwd directed heat fluxes as positive
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_control, only: iday, imonth, isimveg, isimfire
      use inland_comatm
      use inland_comhyd
      use inland_comsno
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_com1d
      use inland_comfire
      use inland_comcrop

      implicit none
!-----------------------------------------------------------------------
! input-output variables
      integer loopi     ! index of little vector in big vector
      integer kpti      ! index of 1st point of little vector 
                        ! in big lpt vector
      integer kptj      ! index of last point of little vector

      integer mcsec     ! current second in day (passed in, 'time' variable)

! local variables
      integer i, k, np, nmtimes
      real*8 rwork, rwork2, rwork3, rwork4, &
             rdepth,    & ! 1/total soil depth over 4 1st layers
             solartot,  & ! total incoming radiation (direct + diffuse, visible + nearIR)
             reflectot, & ! total monthly reflected / monthly incident
             soiltemp,  & ! average soil temp for 4 1st layers
             soilmois,  & ! average soil moisture for 4 1st layers 
             soilice,   & ! average soil ice for 4 1st layers 
             vwc,       & ! total liquid + ice content of 4 1st layers
             awc,       & ! total available water (+ ice) content of 4 1st layer
             snodpth      ! total snow depth
      
      logical lastts ! is this the last timestep of a given day?

!-----------------------------------------------------------------------
! if the first timestep of the month then reset averages
      if (mcsec .eq. 0 .and. iday .eq. 1 ) then
         nmtime(loopi) = 0

         ! TODO reset other vars as well
         amburnfrac(:) = 0
         amburnpft(:,:) = 0
      end if

      lastts = .false.
      if ( mcsec .eq. (86400 - dtime) ) lastts = .true.

! accumulate terms
      nmtime(loopi) = nmtime(loopi) + 1

! working variables
      nmtimes = nmtime(loopi)

! rwork4 for conversion of nitrogen mineralization (moles)
      rwork  = 1. / float(nmtimes)
      rwork2 = float(ndaypm(imonth)) * 86400.
      rwork3 = float(ndaypm(imonth)) * 86400. * 12.e-3
      rwork4 = float(ndaypm(imonth)) * 86400. * 14.e-3
      rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))

!-----------------------------------------------------------------------
! begin global grid
      do i = kpti, kptj

! monthly average temperature
! Different from offline IBIS where average T is from climatology
! TODO: Remove this variabel amts2 - CAGJR
         amts2(i) = ((nmtimes-1) * amts2(i) + ts2(i)) * rwork
      end do
      do 100 i = kpti, kptj

! ---------------------------------------------------------------------
! * * * monthly water budget terms * * *
! ---------------------------------------------------------------------
         amrain(i) = ((nmtimes-1) * amrain(i) + raina(i) * 86400.) * rwork
         amsnow(i) = ((nmtimes-1) * amsnow(i) + snowa(i) * 86400.) * rwork
         amaet(i)  = ((nmtimes-1) * amaet(i) - fvapa(i) * 86400.) * rwork
         amtransu(i) = ((nmtimes-1) * amtransu(i) + gtransu(i) * 86400.) * rwork
         amtransl(i) = ((nmtimes-1) * amtransl(i) + gtransl(i) * 86400.) * rwork
         amsuvap(i) = ((nmtimes-1) * amsuvap(i) + gsuvap(i) * 86400.) * rwork
         aminvap(i) = ((nmtimes-1) * aminvap(i) + ginvap(i) * 86400.) * rwork
         amtrunoff(i) = ((nmtimes-1) * amtrunoff(i) + (grunof(i) + gdrain(i)) *&
                        86400.) * rwork
         amsrunoff(i) = ((nmtimes-1)*amsrunoff(i) + grunof(i)*86400.) *rwork
         amdrainage(i) = ((nmtimes-1)*amdrainage(i) + gdrain(i)*86400.) * rwork
       if(isimagro .gt.0 ) then
         amtrans(i) = ((nmtimes-1) * amtrans(i)  + (gtransl(i) + &
                        gtransu(i)) * 86400.) * rwork
         amtratio(i)  = max(0.0, min(1.0, amtrans(i) / amaet(i))) 

! ---------------------------------------------------------------------
! * * * nitrogen variables * * *
!       converted from kg/m2/s rate to kg/ha/day average
! ---------------------------------------------------------------------
        amtotnleach(i)=((nmtimes-1)*amtotnleach(i)+fout(i,5)/dtime * 1e+04 * 86400.) * rwork
        amno3leach(i)=((nmtimes-1)*amno3leach(i)+nout(i,5)/dtime * 1e+04 * 86400.) * rwork
! ---------------------------------------------------------------------
! * * * monthly atmospheric terms * * *
! ---------------------------------------------------------------------

!        amtemp(i)  = ((nmtimes-1) * amtemp(i)  + ta(i) - 273.16)  * rwork
!        amcloud(i) = ((nmtimes-1) * amcloud(i) + cloud(i) * 100.) * rwork
!        amrh(i)    = ((nmtimes-1) * amrh(i)    + rh(i))           * rwork

      endif !end the isimagro

        amqa(i) = ((nmtimes-1) * amqa(i) + qa(i)) * rwork

! ---------------------------------------------------------------------
! * * * energy budget terms * * *
! ---------------------------------------------------------------------
         solartot = solad(i,1) + solad(i,2) + solai(i,1) + solai(i,2)
         amsolar(i) = ((nmtimes-1) * amsolar(i) + solartot) * rwork
         amirup(i) = ((nmtimes-1) * amirup(i) + firb(i)) * rwork
         amirdown(i) = ((nmtimes-1) * amirdown(i) + fira(i))  * rwork
         amsens(i) = ((nmtimes-1) * amsens(i) - fsena(i)) * rwork
         amlatent(i) = ((nmtimes-1) * amlatent(i) - fvapa(i) * hvap) * rwork

! Total reflected radiation used for albedo calculations
         reflectot = asurd(i,1) * solad(i,1) + asurd(i,2) * solad(i,2) + &
                     asuri(i,1) * solai(i,1) + asuri(i,2) * solai(i,2)

!       amalbedo(i) = ((nmtimes-1) * amalbedo(i) + 
!    >                albedotot/(solartot + epsilon)) * rwork
         amreflect(i) = ((nmtimes-1) * amreflect(i) + reflectot) * rwork

! ---------------------------------------------------------------------
! * * * monthly vegetation parameters * * *
! ---------------------------------------------------------------------
         amlaiu(i) = ((nmtimes-1) * amlaiu(i) + fu(i) * lai(i,2)) * rwork
         amlail(i) = ((nmtimes-1) * amlail(i) + fl(i) * lai(i,1)) * rwork

! ---------------------------------------------------------------------
! * * * monthly soil parameters * * *
! ---------------------------------------------------------------------
         soiltemp = 0.0
         soilmois = 0.0
         soilice  = 0.0
         vwc = 0.0
         awc = 0.0

! averages for first 4 layers of soil (assumed to add to 1 meter depth)
         do 110 k = 1, 4
            soiltemp = soiltemp + tsoi(i,k)  * hsoi(k)
            soilmois = soilmois + wsoi(i,k)  * hsoi(k)
            soilice  = soilice  + wisoi(i,k) * hsoi(k)
            vwc=vwc + (wisoi(i,k) + (1.- wisoi(i,k)) * wsoi(i,k)) * hsoi(k) * &
                poros(i,k)
            awc=awc + max(dble(0.0), (wisoi(i,k) + (1.- wisoi(i,k)) * wsoi(i,k)) -  &
                swilt(i,k)) * hsoi(k) * poros(i,k) * 100.0

110      continue
         soiltemp = soiltemp * rdepth - 273.16
         soilmois = soilmois * rdepth
         soilice  = soilice  * rdepth
         vwc = vwc * rdepth
         awc = awc * rdepth

! monthly average soil parameters:
         amtsoi(i)  = ((nmtimes-1) * amtsoi(i)  + soiltemp) * rwork
         amwsoi(i)  = ((nmtimes-1) * amwsoi(i)  + soilmois) * rwork
         amwisoi(i) = ((nmtimes-1) * amwisoi(i) + soilice)  * rwork
         amvwc(i)   = ((nmtimes-1) * amvwc(i)   + vwc)      * rwork
         amawc(i)   = ((nmtimes-1) * amawc(i)   + awc)      * rwork

! Monthly averages per layer
!
!       do k = 1, nsoilay
!          amtsoil(i,k) = ((nmtimes-1)*amtsoil(i,k)  
!    <                      + tsoi(i,k)) *rwork
!          amwsoil(i,k) = ((nmtimes-1)*amwsoil(i,k)   
!    <                      + wsoi(i,k)) * rwork
!          amwisoil(i,k) = ((nmtimes-1)*amwisoil(i,k)  
!    <                       + wisoi(i,k)) * rwork
!       end do
! ---------------------------------------------------------------------
! * * * snow parameters * * *
! ---------------------------------------------------------------------
         snodpth = hsno(i,1) + hsno(i,2) + hsno(i,3)
         amsnod(i) = ((nmtimes-1) * amsnod(i) + snodpth) * rwork
         amsnof(i) = ((nmtimes-1) * amsnof(i) + fi(i))   * rwork

! ---------------------------------------------------------------------
! * * * fire parameters * * *
! ---------------------------------------------------------------------
if ( isimveg.gt.0 .and. isimfire.eq.2 ) then
         ampbio(i) = ((nmtimes-1) * ampbio(i) + pbio(i)) * rwork
         ampmoi(i) = ((nmtimes-1) * ampmoi(i) + pmoi(i)) * rwork
         ampign(i) = ((nmtimes-1) * ampign(i) + pign(i)) * rwork
         ampfire(i) = ((nmtimes-1) * ampfire(i) + pfire(i)) * rwork
         amsrate(i) = ((nmtimes-1) * amsrate(i) + srate(i)) * rwork
         amabday(i) = ((nmtimes-1) * amabday(i) + abday(i)) * rwork

         ! only update burnfrac at last daily timestep
         !amburnarea(i) = ((nmtimes-1) * amburnarea(i) + burnarea(i)) * rwork
         !amburnpft(i,np) = ((nmtimes-1) * amburnpft(i,np) + burnpft(i,np)) * rwork
         if ( lastts ) then
            amburnfrac(i) = amburnfrac(i) + adburnfrac(i)
            do np = 1, npft
               amburnpft(i,np) = amburnpft(i,np) + adburnpft(i,np)
            end do
         end if

end if
! ---------------------------------------------------------------------
! * * * determine monthly npp * * *
! ---------------------------------------------------------------------
         amnpp(i,1)  = ((nmtimes-1) * amnpp(i,1)  + tnpp(i,1)  * rwork3) * rwork
         amnpp(i,2)  = ((nmtimes-1) * amnpp(i,2)  + tnpp(i,2)  * rwork3) * rwork
         amnpp(i,3)  = ((nmtimes-1) * amnpp(i,3)  + tnpp(i,3)  * rwork3) * rwork
         amnpp(i,4)  = ((nmtimes-1) * amnpp(i,4)  + tnpp(i,4)  * rwork3) * rwork
         amnpp(i,5)  = ((nmtimes-1) * amnpp(i,5)  + tnpp(i,5)  * rwork3) * rwork
         amnpp(i,6)  = ((nmtimes-1) * amnpp(i,6)  + tnpp(i,6)  * rwork3) * rwork
         amnpp(i,7)  = ((nmtimes-1) * amnpp(i,7)  + tnpp(i,7)  * rwork3) * rwork
         amnpp(i,8)  = ((nmtimes-1) * amnpp(i,8)  + tnpp(i,8)  * rwork3) * rwork
         amnpp(i,9)  = ((nmtimes-1) * amnpp(i,9)  + tnpp(i,9)  * rwork3) * rwork
         amnpp(i,10) = ((nmtimes-1) * amnpp(i,10) + tnpp(i,10) * rwork3) * rwork
         amnpp(i,11) = ((nmtimes-1) * amnpp(i,11) + tnpp(i,11) * rwork3) * rwork
         amnpp(i,12) = ((nmtimes-1) * amnpp(i,12) + tnpp(i,12) * rwork3) * rwork

         amnpptot(i) = amnpp(i,1) + amnpp(i,2)  + amnpp(i,3)  + amnpp(i,4) + &
                       amnpp(i,5) + amnpp(i,6)  + amnpp(i,7)  + amnpp(i,8) + &
                       amnpp(i,9) + amnpp(i,10) + amnpp(i,11) + amnpp(i,12)

         if (isimagro .gt. 0) then
             amnpp(i,13)  = ((nmtimes-1) * amnpp(i,13)  + tnpp(i,13)  * rwork3) * rwork
             amnpp(i,14) = ((nmtimes-1) * amnpp(i,14) + tnpp(i,14) * rwork3) * rwork
             amnpp(i,15) = ((nmtimes-1) * amnpp(i,15) + tnpp(i,15) * rwork3) * rwork
             amnpp(i,16) = ((nmtimes-1) * amnpp(i,16) + tnpp(i,16) * rwork3) * rwork
         
             amnpptot(i) = amnpp(i,1) + amnpp(i,2)  + amnpp(i,3)  + amnpp(i,4) +  &
                           amnpp(i,5) + amnpp(i,6)  + amnpp(i,7)  + amnpp(i,8) +  &
                           amnpp(i,9) + amnpp(i,10) + amnpp(i,11) + amnpp(i,12) + &
                           amnpp(i,13)+ amnpp(i,14) + amnpp(i,15)+ amnpp(i,16)
          endif ! end the isimagro

! ---------------------------------------------------------------------
! * * * monthly biogeochemistry parameters * * *
! ---------------------------------------------------------------------
!
! increment monthly total co2 respiration from microbes
! tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
         amco2mic(i) = ((nmtimes-1) * amco2mic(i) + tco2mic(i) * rwork3) * rwork

! increment monthly total co2 respiration from roots
! tco2root is instantaneous value of co2 flux calculated in stats.f
         amco2root(i) = ((nmtimes-1)*amco2root(i) + tco2root(i)*rwork3) * rwork

! calculate average total co2 respiration from soil
         amco2soi(i) = amco2root(i) + amco2mic(i)

!  calculate ratio of root to total co2 respiration
         if (amco2soi(i).gt.0.0) then
            amco2ratio(i) = amco2root(i) / amco2soi(i)
         else
            amco2ratio(i) = -999.99
         endif

!  monthly net ecosystem co2 flux -- npp total minus microbial respiration 
!  the npp total includes losses from root respiration

        if(isimagro .eq. 0) then
            amneetot(i) = ((nmtimes-1) * amneetot(i) + tneetot(i) * rwork3) * rwork
        else
            amneetot(i)  = amnpptot(i) - amco2mic(i)
        endif

! increment monthly total of net nitrogen mineralization
! value for tnmin is calculated in biogeochem.f
         amnmintot(i) = ((nmtimes-1) * amnmintot(i) + &
                        tnmin(i) * rwork4) * rwork
100   continue

      return
end subroutine summonth
