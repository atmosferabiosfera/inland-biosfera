#include "inland_config.h"
#include "inland_compar.h"
! ---------------------------------------------------------------------
subroutine phenocrop (kpti,kptj)
! ---------------------------------------------------------------------
! 
! subroutine which determine the phenological development of each crop
! type including allocation changes of photosynthate to 
! leaf, stem, root, and reproductive organs (pod, seed, grain) 
! leaf area index increase and decline
! subsequent carbon in biomass pools (leaf,stem,root,wood,reproductive)  
!
! Conversion factors for carbon to dry matter are from Penning DeVries
! et al., 1983 and Penning DeVries et al., 1989
! 
! values are fraction of carbon in dry matter
!
! leaves: 0.459
! stem  : 0.494
! root  : 0.467
!
! maize cob and grain : 0.491
! soybean pod and seed: 0.527 
!
! all of the phenology changes are based on the total number of gdd needed
! to change to the next phase - based on fractions of the total gdd typical
! for that region based on the April 1 - Sept 30 window of development  
!
! Phase 1: Planting to leaf emergence
! Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
! Phase 3: Grain fill to physiological maturity and subsequently harvest (LAI decline)  
!
!
! common blocks
!
      use inland_parameters
      use inland_control, only:iyear, iyear0, imonth, iday, jday
      use inland_comatm
      use inland_comsoi
      use inland_comsum
      use inland_comveg
      use inland_comcrop
      use inland_comnitr
      use inland_comsno
      use inland_compft
!
!     real*8 
!	     huigrain(npft), &   ! heat unit index needed to reach vegetative maturity 
!	     huileaf(npft),  &   ! heat unit index needed to attain leaf emergence after planting 
             
!            laicons(npft),  &   ! constant used in lai senescence equation
!            allconsl(npft), &   ! constant used in dynamic allocation equations for leaves (decline in allocation)
!            allconss(npft), &   ! constant used in dynamic allocation equations for stem (decline in allocation) 
!            laimx(npft),    &   ! maximum lai designated for crops from EPIC and EPICphase  models
!            tkill(npft),    &   ! minimum daily temperature threshold used to kill crops
!            mxmat(npft),    &   ! maximum length of growing season to harvest/physiological maturity 
!            mxdgfi(npft),   &   ! maximum number of days needed to reach grain fill stage
!            mxgddgf(npft)       ! maximum gdd past initiation of grain fill to reach phys maturity/harvest
!
!      real*8 laidecl(npoi, npft), &   ! decline in leaf area for crop
!            tauroot(npft),      &   ! time constant for root turnover 
!            arooti(npft),       &   ! initial allocation to crop fine roots 
!            arootf(npft),       &   ! end of growing season final allocation to crop fine roots
!            aleaff(npft),       &   ! end of growing season final allocation to crop leaf area
!            astemf(npft),       &   ! end of growing season final allocation to crop stem biomass
!            fleafi(npft),       &   ! initial fraction of aboveground carbon allocation (stem/leaf) allocated to leaf 
!           fleaf(npft),        &   ! fraction of aboveground carbon allocation (stem/leaf) allocated to leaf 
!           declfact(npft)          ! factor helping to control LAI decline at end of season
!
      real*8 pc,leaftemp,  &
             ddays,        &
             ddfac,        &
             tthreshold,   &
!            xminlai,      &
             dtt,          &
             explf,        &
             tix,          &
             xn,           &
             plag,         &
             abiol,        &
             crmcorn,      &
             crmsgc,sipf3, &
             sipf4,sipf6,  &
             aspecla,ccf5

      integer nplants,     &
              i,j,k,l,n  

! phenology for additional leaf drop - if drought related or temperature related at
! end of growing season     
!
       ddays      = 7.0
       ddfac      = 1.0 / ddays
       tthreshold = 273.16
!
!
       nplants    = 7 
!
      call vernalization(jday)
!
! begin global grid
!      open(900,file='teste_var',status='unknown')
!      write (900,*) "iyear, imonth, iday", iyear, imonth, iday

      do 100 i = lbeg, lend
!
! TODO fix other vars like these, which should be reset every grid point
! !gabriel apagar
! write(*,*) "shape(aleaf)",shape(aleaf)
! write(*,*) "shape(astem)",shape(astem)
! write(*,*) "shape(aroot)",shape(aroot)
! write(*,*) "shape(arepr)",shape(arepr)
! write(*,*) "shape(aerial)",shape(aerial)

!         aleaf(scpft:) = 0.
!         astem(scpft:) = 0.
!         aroot(scpft:) = 0.
!         arepr(scpft:) = 0.
!         aerial(scpft:) = 0.

         if (croptype(i) .ge. scpft) then
            aplantn(i) = 0.0

! only want to look at top 1.5 meters for plant available inorganic n
!
            do 115  k = 1, nsoilay 
               aplantn(i) = aplantn(i) + smsoil(i,k) + smsoln(i,k)
115         continue


!
! Phase 1:
! ==========
!
            huileaf(13)  = lfemerg(13)  * gddmaturity(i,13)
            huileaf(14)  = lfemerg(14)  * gddmaturity(i,14)
            huileaf(15)  = lfemerg(15)  * gddmaturity(i,15)  ! typically between 3-7% in wheat
            if(cropy(i).eq.1) then
               huileaf(16)  = lfemerg(16)  * gddmaturity(i,16)
            else
               huileaf(16)  = lfemerg(16)*(1./3.)* gddmaturity(i,16)
            endif

! Phase 2:
! ==========

            crmcorn      = max(73., min((gddmaturity(i,14)+ 53.683)/13.882,135.))

            huigrain(14) = -0.002  * (crmcorn - 73.) + grnfill(14)
            huigrain(14) = min(max(huigrain(14),grnfill(14) - 0.1), grnfill(14)) 
            huigrain(14) = huigrain(14)   * gddmaturity(i,14)  ! from Cabelguenne et al. 1999

            crmsgc       = max(73., min((gddmaturity(i,16)+ 53.683)/13.882,135.))

            huigrain(16) = -0.002  * (crmsgc - 73.) + grnfill(16)
            huigrain(16) = min(max(huigrain(16),grnfill(16) - 0.1), grnfill(16)) 
            huigrain(16) = huigrain(16)   * gddmaturity(i,16)  

            huigrain(13) = grnfill(13)    * gddmaturity(i,13)  ! from Cabelguenne et al. 1999

! TODO FIXME here and in other places also...
            if ( iwheattype .ne. 0 ) then
               huigrain(15) = grnwht(iwheattype) * gddmaturity(i,15)
               fleafi(15)   = fleafiwht(iwheattype)  ! wheat          
               mxgddgf(15)  = mgddgf(iwheattype)  ! wheat
               mxmat(15)    = mxmatwht(iwheattype)   ! wheat
            else
               huigrain(15) = 0.0
               fleafi(15)   = 0.0
               mxgddgf(15)  = 0.0
               mxmat(15)    = 0.0
            endif

            do 80 j = scpft, ecpft 
               cropout(i,j,31) = 0 !aleaf(j)
               cropout(i,j,32) = 0 !aroot(j)
               cropout(i,j,33) = 0 !arepr(j)
               cropout(i,j,34) = 0 !astem(j)


               if (croplive(i,j) .eq. 1.0) then
                  grnfraccrop(i,j) = 1.0

                  fleaf(j) = fleafi(j) * (exp(-bfact(j)) - exp(-bfact(j) * &
                             gddplant(i,j)/huigrain(j))) / (exp(-bfact(j))-1)
!
! calculate accumulated growing degree days since planting (gddplant) 
! determine if growing degree days calculated from top layer soil temperature
! are enough for leaf emergence to occur 

                  hui(i,j) = gddplant(i,j)
   
                  if(j.eq.16) then
                     leafout(i,j) = gddplant(i,j)
                  else
                     leafout(i,j) = gddtsoi(i,j) 
                  endif
                  laidecl(i,j) = 0.0
                  idpp(i,j) = idpp(i,j) + 1
!
! crop phenology from leaf emergence to start of leaf decline   
! determine allocation coefficients to help calculate increase in lai  
! and fine root biomass
!
                  if (leafout(i,j) .ge. huileaf(j) .and. &
                      hui(i,j) .lt. huigrain(j).and.j.ne.16) then

! Phase 1 completed:
! ==================
                     awood(i,j) = 0.0
!           
! check to see if lai is at maximum allowable 
!
                  aroot(i,j) = min(1.0, (arooti(j) - (arooti(j) - arootf(j)) * &
                               min(1.0,hui(i,j)/gddmaturity(i,j))))

                  aroot(i,j) = max(0.0, aroot(i,j))

                     if (peaklai(i,j) .eq. 1) then 
                        aleaf(i,j) = 0.0
                        arepr(j) = 0.0
                        astem(j) = 1.0 - aroot(i,j) - aleaf(i,j) - arepr(j)
                     else
                        aleaf(i,j) = max(0.0,(1.0 - aroot(i,j)) * fleaf(j))
                        astem(j) = 1.0 - aroot(i,j) - aleaf(i,j)
                        arepr(j) = 0.0
                     endif

                     if (peaklai(i,j) .eq. 0) then
                        tlai(i,j) = plai(i,j) + (specla(i,j) * aleaf(i,j) * &
                                    max(0.0, adnpp(i,j)))

                        if (tlai(i,j) .ge. laimx(j)) then
                           aleaf(i,j) = min(1.0 - aroot(i,j), laimx(j) - &
                                      plai(i,j)) / (specla(i,j) * adnpp(i,j))
                           aleaf(i,j) = max(0.0, aleaf(i,j))
                           astem(j) = 1.0 - aroot(i,j) - aleaf(i,j)

!
! CJK other possible source of over allocation
! above astem could be set to 0.0
!
                           peaklai(i,j) = 1
                        endif
!
! original phenology - apply to soybeans
! if maize - employ leaf expansion routine based on a modified
! ceres gdd approach
!
                        plai(i,j) = plai(i,j) + (specla(i,j) * aleaf(i,j) * &
                                    max(0.0, adnpp(i,j)))

                     endif !of peak lai

! hold ending allocation values to stem and leaf for use by equations after shift to
! reproductive phenology stage begins
!
                     astemi(j) = astem(j)
                     aleafi(j) = aleaf(i,j)

!
! shift allocation either when enough gdd are accumulated or maximum number
! of days has elapsed since planting
!
                  else if (hui(i,j) .ge. huigrain(j) .and. j .ne. 16 &
                          .and. croplive(i,j) .eq. 1.) then   
                     dpgf(i,j) = max(0.0, gddplant(i,j) - huigrain(j))
!
! keep day of year that grain fill begins
!
                     grainday(i,j) = min(grainday(i,j), real(jday))
!
! Phase 2 completed:
!
                     awood(i,j) = 0.0
                     aroot(i,j) = min(1.0, (arooti(j) - (arooti(j) - arootf(j))* &
                                min(1.0, hui(i,j)/gddmaturity(i,j)))) 
                     aroot(i,j) = max(0.0, aroot(i,j))


                     if (thrlai(i,j) .lt. 0.0001) thrlai(i,j) = plai(i,j)
!
! lai plateau -- reached threshold -- according to heat unit index  (accumulated GDD logic)        
! set lai and hui threshold  (to be used in leaf phenology (senescence) calculations)
! shift aboveground allocation to grain/fruit
!
                     templai(i,j)   = plai(i,j)

! 
! lai decline based thermal time accumulation past grain fill 
! add check in to make sure huigrain is <= gddplant in this equation.  If it
! isn't because the days past planting is used to shift grain allocation, set
! huigrain value to the gddplant value when the shift took place so LAI starts
! the proper decline.
! 
                     plai(i,j) = max(thrlai(i,j) * (1.0-min(max((gddplant(i,j)-     &
                                 huigrain(j)),0.0)/(0.55*gddmaturity(i,j)), 1.0) ** &
                                 laicons(j)), xminlai)

! calculate decrease in lai for purpose of updating aboveground biomass pools
!
                     laidecl(i,j) = max(0.0, templai(i,j) - plai(i,j))

!
                     if (astemi(j) .gt. astemf(j)) then
                        astem(j) = max(astem(j) * (1.0-min((hui(i,j)-         &
                                huigrain(j))/((gddmaturity(i,j)*declfact(j))- &
                                huigrain(j)),1.0)**allconss(j)),astemf(j)) 
                        astem(j) = max(0.0,astem(j))
                     endif

! CJK 9-23-04
! ***********
write(*,*) "===================================================================PASSEI CANA"
                     if (aleafi(j) .gt. aleaff(j)) then
                        aleaf(i,j) = max(aleaf(i,j) * (1.0-min((hui(i,j)-     &
                            huigrain(j))/((gddmaturity(i,j)*declfact(j))- &
                            huigrain(j)),1.0)**allconsl(j)),aleaff(j)) 
                        aleaf(i,j) = max(0.0,aleaf(i,j)) 
                     endif
                        arepr(j) = 1.0 - aroot(i,j) - astem(j) - aleaf(i,j) - awood(i,j)
                  endif

!                  if((iyear.eq.2009.and.jday.gt.300).or.(iyear.eq.2010.and.jday.lt.150)) then
                  if((imetyear .ne. 9999)) then
                        write(42,422) iyear,jday, idpp(i,j),hui(i,j), aroot(i,j),aleaf(i,j), &
                                      astem(j),arepr(j),plai(i,j)  
                  endif
422  format(2(i4,1x),9(f7.2,1x))

                  if (leafout(i,j).lt.huileaf(j).and.j.eq.16)  then
                     gddemerg(i) = 0.0
                     awood(i,j) = 0.0
                     aroot(i,j) = 0.0
                     aerial(j) = 0.0
                     aleaf(i,j) = 0.0
                     astem(j) = 0.0
                     arepr(j) = 0.0
                     rm(i) = 0.0 

                  else if (leafout(i,j).ge.huileaf(j).and.j.eq.16)  then
! Phase 1 completed:
! ==================
!gabriel apagar
write(*,*) "===================================================================PASSEI CANA"
                     if(gddemerg(i).eq.0) gddemerg(i) = gddplant(i,j)

                        awood(i,j) = 0.0
                        aroot(i,j) = 0.00001
                        aerial(j) = 0.00001
                        aleaf(i,j) = 0.00001
                        astem(j) = 0.00001
                        arepr(j) = 0.00001
                        rm(i)= min(100.0, 100 * (gddplant(i,j)-gddemerg(i)) / (gddmaturity(i,j)-gddemerg(i)) )
                        rm(i)=max(0.000000000000001,rm(i))

                        if(cropy(i).eq.1) then
                           aerial(j) = ( 1-arootf(j) ) * min(1.0, ( 1-exp(-rootd*0.2*rm(i)) ) )
                        else
                           aerial(j) = ( 1-arootf(j) ) * min(1.0, ( 1-exp(-rootd*    rm(i)) ) ) 
                        endif
                        aerial(j)=max(0.000000000000001,aerial(j))
                        
                        aerial(j) = ( 1-arootf(j) ) * min(1.0, ( 1-exp(-rootd*rm(i)) ) ) 

                        aroot(i,j) = 1.-aerial(j)
                        af1 = max(0.0, rm(i)*sf1 - sf1*ipf1)
                        af2 = max(0.0, 1.0 - ( exp(ecf2*ipf2)/exp(ecf2*rm(i)) ) )

                        astem(j) = aerial(j) * min(1.0, max(af1,af2))
                        astem(j) = min( max(0.1,aerial(j)-aleaff(j)) ,astem(j) )

                        aleaf(i,j) = aerial(j)-astem(j)
                        astem(j) = astem(j) + min( aleaf(i,j), ldf*aleaf(i,j)* &
                                   min(1.0,exp(tmld*ecf7)/exp((td(i)-273.16)*ecf7) ) )


                        astem(j) = min( max(0.0,aerial(j)-aleaff(j)) ,astem(j) )
                        aleaf(i,j) = aerial(j)-astem(j)
                        sipf3=ipf1+(100.-ipf1)*(ipf3/100.)
                        af3 = max(0.0, rm(i)*sf3 - sf3*sipf3)
                        sipf4=ipf1+(100.-ipf1)*(ipf4/100.)
                        af4 = max(0.0, 1.0- ( exp(ecf4*sipf4)/exp(ecf4*rm(i)) ) )
                        arepr(j) = astem(j) * min(1.0, max(af3,af4))
                        arepr(j) = min(aerial(j)-aleaf(i,j), arepr(j) )
                        astem(j)= astem(j)- arepr(j)
                        af5= min(1., max(0., 1. -( exp((td(i)-273.16)*ecf5)/exp(tf5*ecf5) ) ) ) + &
                             min(0., min(0., ( exp(tf5*ecf5)/exp((td(i)-273.16)*ecf5) ) -1 ) )
                        sipf6=sipf4+(100.-sipf4)*(ipf6/100.)
                        af6 = max(0.0, 1.0- ( exp(ecf6*sipf6)/exp(ecf6*rm(i)) ) )
                        ccf5=arepr(j)
                        arepr(j) = arepr(j)+ astem(j)*wf5*af5*af6
                        arepr(j) = min(aerial(j)-aleaf(i,j), arepr(j) )
                        astem(j)= astem(j) - (arepr(j) - ccf5)       
                        leaftemp=aleaf(i,j)
                        tlai(i,j) = plai(i,j) + (specla(i,j) * aleaf(i,j) * &
                                    max(0.0, adnpp(i,j)))

                        if (tlai(i,j) .ge. laimx(j)) then
                           aleaf(i,j) = min( aleaf(i,j),(laimx(j)-plai(i,j)) / &
                                      (specla(i,j) * adnpp(i,j)))
                           aleaf(i,j) = max(0.0, aleaf(i,j))
                           aroot(i,j) = aroot(i,j) + ( (leaftemp-aleaf(i,j)) * aroot(i,j)/(astem(j)+arepr(j)+aroot(i,j)) )
                           astem(j) = astem(j) + ( (leaftemp-aleaf(i,j)) * astem(j)/(astem(j)+arepr(j)+aroot(i,j)) )
                           arepr(j) = arepr(j) + ( (leaftemp-aleaf(i,j)) * arepr(j)/(astem(j)+arepr(j)+aroot(i,j)) )
                        endif  
                     aroot(i,j) = max(0.0, aroot(i,j))
                     aleaf(i,j) = max(0.0, aleaf(i,j))
                     astem(j) = max(0.0, astem(j))
                     arepr(j) = max(0.0, arepr(j))

                     if (arepr(j) .gt. 0.001 .and. grainday(i,j) .gt. 999) grainday(i,j) = min(grainday(i,j), &
                        real(jday))
                        templai(i,j) = (cbiol(i,j)*specla(i,j))
                        plai(i,j) = (cbiol(i,j)*specla(i,j)) - (cbiol(i,j)*specla(i,j))* (  (1./tauleaf(j)) )
                 
                     if(td(i) .le. 278.16 .and. td(i) .ge. 268.16) then
                        plai(i,j) = plai(i,j)* max(0.4,min(1.,0.5+((td(i)-268.16)/20.) )) 
                     endif
                     plai(i,j) = max(0.01,plai(i,j))
                     laidecl(i,j)   = max(0.0, templai(i,j) - plai(i,j))

                  endif

! keep track of total biomass production for the entire year, and the
! aboveground value to calculate harvest index

                  aybprod(i,j) = aybprod(i,j)                    &
                               + aleaf(i,j) * max(0.0,adnpp(i,j))  &
                               + astem(j) * max(0.0,adnpp(i,j))  &
                               + aroot(i,j) * max(0.0,adnpp(i,j))  &
                               + awood(i,j) * max(0.0,adnpp(i,j))  &
                               + arepr(j) * max(0.0,adnpp(i,j))
!-above
                  ayabprod(i,j) = ayabprod(i,j)                  & 
                                + aleaf(i,j) * max(0.0,adnpp(i,j)) &
                                + astem(j) * max(0.0,adnpp(i,j)) &
                                + arepr(j) * max(0.0,adnpp(i,j)) &
                                + awood(i,j) * max(0.0,adnpp(i,j))

                  cropout(i,j,31) = aleaf(i,j)
                  cropout(i,j,33) = arepr(j)

!
! keep track of annual total root production carbon
!
                  ayrprod(i,j) = ayrprod(i,j) + aroot(i,j) * max(0.0,adnpp(i,j))
!
! keep track of total carbon allocated to leaves for litterfall calculation
!
                  aylprod(i,j) = aylprod(i,j) + aleaf(i,j) * max (0.0, adnpp(i,j))
!

                  cbiol(i,j) = cbiol(i,j) + aleaf(i,j) * max (0.0, adnpp(i,j)) - &
                               (laidecl(i,j)/specla(i,j))
                  cbiog(i,j) = cbiog(i,j) + arepr(j) * max (0.0, adnpp(i,j))
                  cbios(i,j) = cbios(i,j) + astem(j) * max (0.0, adnpp(i,j))

                  fallrsgc(i,3) = cbior(i,j) + aroot(i,j) * max(0.0,adnpp(i,j))
                  cbior(i,j) = cbior(i,j) * exp(-1.0/tauroot(j)) +             &
                                aroot(i,j) * tauroot(j) * max(0.0,adnpp(i,j)) *  &
                                (1.0 - exp(-1.0/tauroot(j)))
                  fallrsgc(i,3) = fallrsgc(i,3) - cbior(i,j)
                  falllsgc(i,3) = laidecl(i,j)/specla(i,j)
                  cbiow(i,j) = 0.0 
                  cbior(i,j) = max(0.0, cbior(i,j))
                  cbiol(i,j) = max(0.0, cbiol(i,j)) 
                  cbios(i,j) = max(0.0, cbios(i,j))
                  cbiog(i,j) = max(0.0, cbiog(i,j))

!
! update vegetation's physical characteristics
!
                  plai(i,j)    = cbiol(i,j) * specla(i,j) 

                  if(j.eq.16.and.rm(i).gt.1.) then
                     if(cbiol(i,j).eq.0)then
                        cbiol(i,j)=0.000001
                     endif
                     plai(i,j) = cbiol(i,j)*specla(i,j) + (aylprod(i,j)- cbiol(i,j) ) *specla(i,j)* 0.15 
                     grnfraccrop(i,j)= cbiol(i,j)*specla(i,j)/plai(i,j) 
                  endif
                  biomass(i,j) = cbiol(i,j) + cbiog(i,j) + cbior(i,j)+  &
                                 cbios(i,j) + cbiow(i,j)
!
! keep track of aboveground annual npp 
!
                  ayanpp(i,j) = (aleaf(i,j) + arepr(j) + astem(j) + &
                                awood(i,j)) * adnpp(i,j) + ayanpp(i,j)

       if(imetyear .ne. 9999) then        
!sant - convertendo de kg.C/m2 para  tonelada MS/ha
          if(npoi .eq. 1)then
            write(222,42)iyear,jday,idpp(i,j),plai(i,j),plai(i,j)*grnfraccrop(i,j),             &
                         cbiol(i,j)*(1/0.45)*10.0,(aylprod(i,j)-cbiol(i,j))*(1/0.45)*10.0,      &
                         (aybprod(i,j)-ayrprod(i,j)-aylprod(i,j)-cbiog(i,j))*(1/0.45)*10.0,     &
                         cbiog(i,j)*(1/0.45)*10.0,(aybprod(i,j)-ayrprod(i,j))*(1/0.45)*10.0,    &
                         ayrprod(i,j)*(1/0.45)*10.0
         
!sant- usando para o perfil de raiz	qgalho=0
          endif
 42     format (i4,1x,i3,1x,9(f6.2,1x)) 
	endif

!---------------------------------------------------------------------------------
! check for climatic and phenological limits on maturity, growth, and harvest date
!
! check to see if minimum temperature has fallen below freeze
! kill threshold for 3 consecutive days and if lai is above a minimum, plant will
! be damaged/killed.  This function is more for spring freeze events
! or for early fall freeze events
!
! currently simulates too many grid cells that are killed by
! freezing temperatures 
!
! spring wheat is affected by this, winter wheat kill function
! is determined in crops.f - is a more elaborate function of
! cold hardening of the plant
!
                  if (tmin(i) .le. tkill(j)) then
                     ccdays(i,j) = ccdays(i,j) + 1
                  else
                     ccdays(i,j) = 0
                  endif
!
                  if (ccdays(i,j) .ge. 1 .and.                       &
                     hui(i,j) .ge. 0.6*gddmaturity(i,j) .and.        &
                     croplive(i,j) .eq. 1 .and.                      &
                     ((j .eq. 13 .or. j .eq. 14 .or. j .eq. 16) .or. &
                     j .eq. 15 .and. iwheattype .eq. 1)) then
                        croplive(i,j) = 0.0
                        harvdate(i,j) = jday
                  endif

                  if (j .eq. 16) then   
                     if (cropy(i) .eq. 1) then
                        if (((hui(i,j) .ge. gddmaturity(i,j)) .and. (idpp(i,j) .ge.        &
                           mxmat(j)-60) .and. (iday .ge. pdmin(j) .and. imonth .ge.        &
                           (mod(pmmin(16)+(mxmat(16)/30)-3,12.0)+1))) .or. ((idpp(i,j) .ge. &
                           mxmat(j)-30) .and. (iday .eq. pdmin(j) .and. imonth .eq.        &
                           (mod(pmmin(16)+(mxmat(16)/30)-1,12.0) + 1)))) then
                              croplive(i,j) = 0.0
                              grnfraccrop(i,j) = 0.0           
                              if (harvdate(i,j) .eq. 999) harvdate(i,j) = jday 
                              plai(i,j) = 0.25      
                        endif
                     else
                        if (((hui(i,j) .ge. gddmaturity(i,j)) .and. (idpp(i,j) .ge. 335))  &
                            .or. idpp(i,j) .ge. 395 .or. ((idpp(i,j) .ge. 335) .and.       &
                           (iday .eq. pdmin(j) .and. imonth .eq. (mod(pmmin(16) +          &
                           (mxmat(16)/30)-1,12.0)+1)))) then             !maximum expected harvest date
                              croplive(i,j) = 0.0
                              grnfraccrop(i,j) = 0.0           
                              if (harvdate(i,j) .eq. 999) harvdate(i,j) = jday 
                              plai(i,j) = 0.25       
                        endif
                     endif
                  else
                     if (hui(i,j) .ge. gddmaturity(i,j) .or. &
                     idpp(i,j) .ge. mxmat(j)) then
                        croplive(i,j) = 0.0
                        grnfraccrop(i,j) = 0.0      
                        if (harvdate(i,j) .eq. 999) harvdate(i,j) = jday 
                           plai(i,j) = 0.25
                     endif
                  endif
               endif
80          continue
         endif

100   continue

        call cropupdate(jday)
!
        call cropresidue(jday) 
!
        if (iday .eq. 31 .and. imonth .eq. 12) then 
          call cropoutput(jday) 
        endif 
!
      return
!
 end subroutine phenocrop
!
