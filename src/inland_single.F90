#include "inland_config.h"

#ifndef SINGLE_POINT_MODEL
#error "This subroutine should ONLY be compiled for 0D INLAND model option."
#endif

subroutine single(linenum, test)

      use inland_comatm, only: fira, solai, asuri, solad, asurd
      use inland_com1d, only: firb, fsena, fvapa, abupi, abupd
      use inland_comveg, only: topparl, topparu, chl, chs, chu, tu, ts, tl, &
                               sai, lai, tco2mic, tco2root, fl, fu, wliql,  &
			       wliqs, wliqu, tnpptot, tgpptot, totlit, totlail, totlaiu, &
			       totcsoi, cbioltot, cbiortot, cbiowtot, cbior, cbiow, cbiol, &
			       clitll, clitlm, clitrm, clitws, clitrs, clitrl, clitwm, falll, &
			       tneetot, clitls, clitwl
      use inland_control, only: nan
      use inland_comsoi, only: tg, hsoi, poros, swilt, soihfl, wsoi, upsoil, tsoi, &
                               rhosoi, csoi
      use inland_parameters, only: hvap, nsoilay, rhow
      use inland_comhyd, only: gsuvap, grunof, gdrain, ginvap, gtrans
      use inland_comsum, only: wtot, ayanpptot
      use inland_comforc, only: xprec, xcld, xlati, xlin, xsin, xta, xqa,    &
                                xua, dimforc
		
      implicit none
	       
      real *8 rn,          & ! net radiation flux (short and long waves)
              pari,        & ! incoming PAR
              apar,        & ! APAR
              paro,        & ! outgoing PAR
	      swnet,       & ! net shortware radiation
              lwnet,       & ! net longwave radiation
              solartot,    & ! total incoming radiation (dir+diff,vis+nIR)
              albedo,      & ! surface albedo
              vegt,        & ! vegetation canopy temperature
              avgsurft,    & ! average surface temperature
              reflsw,      & ! reflected solar radiation
              canopint,    & ! total canopy water storage
              rootmoist,   & ! root zone soil moisture
              soilwet,     & ! total soil wetness
              surfalb,     & ! surface albedo
              totsoilcarb, & ! total soil carbon (kg.C.m^-2)
              autoresp,    & ! autotrophic respiration (kg.m^-2.s^-2)
              surfheat,    & ! surface heat storage (J.m^-2)
              soilmoistsum,& ! soil moisture in all layers (kg.h2o.m^-2)
	      delsurfstor, & ! change in surface heat storage (J.m^-2)
              delwtot,     & ! change in total water stored (kg.h2o.m^-2.s^-1)
	      delsoilmoist,& ! change in soil moisture (kg.h2o.m^-2)
              delsurfheat, & ! change in surface heat storage (J.m^-2)
              surfheatl,   & ! last surface heat storage (J.m^-2)
              soilmoist(nsoilay), & ! soil moisture in which layer (kg.h2o.m^-2)
              soilmoistl,  & ! last soil moisture (kg.h2o.m^-2)
              wtotl,       & ! total wather stored in soil, vegetation and snow
              fpar,        & ! absorbed fraction of PAR
              agb,         & ! above ground biomass
	      totlai,      & ! total ecosystem lai
	      totlivbio,   & ! total living biomass (kg.C.m^-2)
              delcanheat     ! change incanopy heat storage (J.m^-2)
	      
      integer i, linenum     ! line number
      
      real *8 test	     ! test on timestep

! FIXME ET: this should not be inside the main loop, but in a function or subroutine...
! energy budget of the surface


                  rn = solad(1,1) * (1-asurd(1,1)) + &
                       solad(1,2) * (1-asurd(1,2)) + &
                       solai(1,1) * (1-asuri(1,1)) + &
                       solai(1,2) * (1-asuri(1,2)) + fira(1) - firb(1)

                  pari = (solad(1,1) + solai(1,1)) * 4.59
                  apar = (topparl(1) + topparu(1)) + 4.59
                  paro = (solad(1,1) * asurd(1,1) + solai(1,1)*asuri(1,1))*4.59

                  swnet = solad(1,1) * (1-asurd(1,1)) + &
                          solad(1,2) * (1-asurd(1,2)) + &
                          solai(1,1) * (1-asuri(1,1)) + &
                          solai(1,2) * (1-asuri(1,2))

                  lwnet = (fira(1) - firb(1))

                  reflsw = solad(1,1)*asurd(1,1) + solad(1,2)*asurd(1,2) + &
                           solai(1,1)*asuri(1,1) + solai(1,2)*asuri(1,2)

                  solartot = solad(1,1) + solad(1,2) + solai(1,1) + solai(1,2)
		  
! Imbuzeiro: Created to replace the NaN values to zero, because
!            in time whitout radiation the inland write NaN and I
!            excluded the NaN to used the output file in my
!            netcdf program (.txt to .nc)
! In accordance to Hewlley Imbuzeiro, when solartot is 0, albedo must be an 
! 'expected' NaN (or -999.99 defined by the 'nan' variable on inland_control.F90) - fzm

                  if (solartot.ne.0) then
                     albedo = reflsw / solartot
                  else
                     albedo = nan
                  endif

                  vegt = ((chl*tl(1)) + (chs*ts(1)) + (chu*tu(1))) / &
                         (chs+chl+chu)

                  avgsurft = (tu(1) + tl(1) + tg(1)) / 3
		  		  
! Surfalb is -999.99 (expected NaN) when solad+solai = 0 (no solar incidence)
                  if ((solad(1,1)+solai(1,1)).ne.0.) then
                     surfalb = ((solad(1,1)*(1-asurd(1,1))) + (solai(1,1) * &
                               (1-asuri(1,1)))) / (solad(1,1) + solai(1,1))
                  else
                     surfalb = nan
                  endif

                  soilwet = ((wsoi(1,1) - swilt(1,1) * hsoi(1) * poros(1,1)) + &
                            (wsoi(1,2) - swilt(1,2) * hsoi(2) * poros(1,2)) +  &
                            (wsoi(1,3) - swilt(1,3) * hsoi(3) * poros(1,3)) +  &
                            (wsoi(1,4) - swilt(1,4) * hsoi(4) * poros(1,4)) +  &
                            (wsoi(1,5) - swilt(1,5) * hsoi(5) * poros(1,5)) +  &
                            (wsoi(1,6) - swilt(1,6) * hsoi(6) * poros(1,6))) / &
                            ((1-swilt(1,1) * poros(1,1) * hsoi(1)) +           &
                            (1-swilt(1,2) * poros(1,2) * hsoi(2)) +            &
                            (1-swilt(1,3) * poros(1,3) * hsoi(3)) +            &
                            (1-swilt(1,4) * poros(1,4) * hsoi(4)) +            &
                            (1-swilt(1,5) * poros(1,5) * hsoi(5)) +            &
                            (1-swilt(1,6) * poros(1,6) * hsoi(6)))

                  canopint = (wliqu(1) * fu(1) * 2 * lai(1,2) + &
                              wliqs(1) * fu(1) * 2 * sai(1,2) + &
                              wliql(1) * fl(1) * 2 * (lai(1,1) + sai(1,1)))

                  rootmoist = (upsoil(1,1) + gsuvap(1))
                  totsoilcarb = tco2root(1) + tco2mic(1)
                  delcanheat = (swnet + lwnet - fvapa(1)*hvap - fsena(1) + &
                                soihfl(1)) * 3600
				
! Imbuzeiro: Correct the soilmoisture equation
! Turned the code into a do-loop for simplicity - fzm
                  soilmoistsum=0.
                  do 260 i=1,nsoilay
                     soilmoist(i) = (wsoi(1,i) * poros(1,i) * rhow * hsoi(i))
                     soilmoistsum = soilmoistsum + soilmoist(i)
260               continue

                  delsoilmoist = soilmoistsum - soilmoistl
                  soilmoistl = soilmoistsum

                  surfheat=0.
                  do 270 i=1,nsoilay
                     surfheat=surfheat+(csoi(1,i)*rhosoi(1,i)*tsoi(1,i)*hsoi(i))
270               continue

                  delsurfheat = surfheat - surfheatl
                  surfheatl = surfheat
                  
                  delwtot = wtot(1) - wtotl
                  wtotl = wtot(1)

                  autoresp = tgpptot(1) - tnpptot(1)
                  totsoilcarb = totcsoi(1) + totlit(1)

! Imbuzeiro: Created to replace the NaN values to zero, because
!            in time whitout radiation the inland write NaN and I
!            excluded the NaN to used the output file in my
!            netcdf program (.txt to .nc)
! In accordance to Hewlley Imbuzeiro, when solartot is 0, albedo must be an 
! 'expected' NaN (or -999.99 defined by the 'nan' variable on inland_control.F90) - fzm
                  if (solad(1,1).ne.0) then
                     fpar = (solad(1,1)*abupd(1) + solai(1,1)*abupi(1)) / &
                            (solad(1,1) + solai(1,1))
                  else
                     fpar = nan
                  endif

! FIXME: make a DO loop for totlivbio, agb, cbio?tot all at once -fzm
                  totlivbio = cbiol(1,1) + cbior(1,1) + cbiow(1,1) + &
                              cbiol(1,2) + cbior(1,2) + cbiow(1,2) + &
                              cbiol(1,3) + cbior(1,3) + cbiow(1,3) + &
                              cbiol(1,4) + cbior(1,4) + cbiow(1,4) + &
                              cbiol(1,5) + cbior(1,5) + cbiow(1,5) + &
                              cbiol(1,6) + cbior(1,6) + cbiow(1,6) + &
                              cbiol(1,7) + cbior(1,7) + cbiow(1,7) + &
                              cbiol(1,8) + cbior(1,8) + cbiow(1,8) + &
                              cbiol(1,9) + cbior(1,9) + cbiow(1,9) + &
                              cbiol(1,10) + cbior(1,10) + cbiow(1,10) + &
                              cbiol(1,11) + cbior(1,11) + cbiow(1,11) + &
                              cbiol(1,12) + cbior(1,12) + cbiow(1,12)

                  agb = cbiol(1,1) + cbiow(1,1) + cbiol(1,2) + cbiow(1,2) + &
                        cbiol(1,3) + cbiow(1,3) + cbiol(1,4) + cbiow(1,4) + &
                        cbiol(1,5) + cbiow(1,5) + cbiol(1,6) + cbiow(1,6) + &
                        cbiol(1,7) + cbiow(1,7) + cbiol(1,8) + cbiow(1,8) + &
                        cbiol(1,9) + cbiow(1,9) + cbiol(1,10) + cbiow(1,10) + &
                        cbiol(1,11) + cbiow(1,11) + cbiol(1,12) + cbiow(1,12)

                  cbioltot(1) = cbiol(1,1) + cbiol(1,2) + cbiol(1,3) + &
                                cbiol(1,4) + cbiol(1,5) + cbiol(1,6) + &
                                cbiol(1,7) + cbiol(1,8) + cbiol(1,9) + &
                                cbiol(1,10) + cbiol(1,11) + cbiol(1,12)

                  cbiortot(1) = cbior(1,1) + cbior(1,2) + cbior(1,3) + &
                                cbior(1,4) + cbior(1,5) + cbior(1,6) + &
                                cbior(1,7) + cbior(1,8) + cbior(1,9) + &
                                cbior(1,10) + cbior(1,11) + cbior(1,12)

                  cbiowtot(1) = cbiow(1,1) + cbiow(1,2) + cbiow(1,3) + &
                                cbiow(1,4) + cbiow(1,5) + cbiow(1,6) + &
                                cbiow(1,7) + cbiow(1,8) + cbiow(1,9) + &
                                cbiow(1,10) + cbiow(1,11) + cbiow(1,12)

                  totlai = totlaiu(1) + totlail(1)

! Imbuzeiro: For the LBA-DMIP protocol I did a spinup but I only wanted
!            write the last 3 years. This information changes for other
!            sites and user specifications
                     linenum = linenum + 1
                     write(40,9200) linenum, swnet, lwnet, -fvapa(1)*hvap,     &
                                    -fsena(1), soihfl(1), delcanheat, nan,     &
                                    -fvapa(1), grunof(1), nan, gdrain(1),      &
                                    gdrain(1)+grunof(1), delsoilmoist, nan,    &
                                    nan, nan, vegT,tsoi(1,1), nan, albedo, nan,&
                                    fpar,(soilmoist(i),i=1,nsoilay),           &
                                    (tsoi(1,i),i=1,nsoilay), soilwet,ginvap(1),&
                                    gtrans(1), gsuvap(1), nan,nan, canopint,   &
                                    tgpptot(1)*0.012,tnpptot(1)*0.012,         &
                                    -tneetot(1)*0.012,autoresp*0.012,          &
                                    tco2mic(1)*0.012,ayanpptot(1), nan,        &
                                    falll(1),cbioltot(1), cbiortot(1),         &
                                    cbiowtot(1),clitll(1), clitlm(1),          &
                                    clitls(1),clitrl(1), clitrm(1), clitrs(1), &
                                    clitwl(1), clitwm(1), clitws(1), agb,      &
                                    totlivbio, cbiowtot(1), totlai,            &
                                    xta(1,linenum)+273.15, xqa(1,linenum),     &
                                    xua(1,linenum), xprec(1,linenum)/1800, nan,&
                                    xsin(1,linenum), xlin(1,linenum), nan
!                  endif                       
9200              format(i10,x,74(f25.16,x))

                  test = test - 1

! end of half-hourly loop

      return
end subroutine single
