#include "inland_config.h"

subroutine build_file


#ifdef SINGLE_POINT_MODEL

! Initialize the output file from outside the time loops -fzm
! Write half-hourly output
! Imbuzeiro: Added the variable names in the output file
      open(40,file='output/single_point-output.csv',status='unknown')
      write(40,9199) 'time','Swnet','Lwnet','Qle','Qh','Qg','DelCanHeat',   &
                     'DelSurfHeat','Evap','Qs','Qrec','Qsb','Qt',           &
                     'DelSoilMoist','DelSurfStor','DelIntercept',           &
                     'WaterTableD','VegT','BaresoilT','AvgSurfT','Albedo',  &
                     'SurfStor','fPAR','SoilMoist1','SoilMoist2',           &
                     'SoilMoist3','SoilMoist4','SoilMoist5','SoilMoist6',   &
                     'SoilTemp1','SoilTemp2','SoilTemp3','SoilTemp4',       &
                     'SoilTemp5','SoilTemp6','SoilWet','Ecanop','Tveg',     &
                     'Esoil','Ewater','RootMoist','CanopInt','GPP',         &
                     'NPP','NEE','AutoResp','HeteroResp','ANPP','BgResp',   &
                     'LitterFall','CarbPools1','CarbPools2','CarbPools3',   &
                     'CarbPools4','CarbPools5','CarbPools6','CarbPools7',   &
                     'CarbPools8','CarbPools9','CarbPools10','CarbPools11', &
                     'CarbPools12','CarbPools13','TotLivBiom','AbvGrndWood',&
                     'LAI','Tair','Qair','Wind','Rainf','Psurf','Swdown',   &
                     'Lwdown','CO2air'
9199  format (74(a,x))

#else /* SINGLE_POINT_MODEL */

          open (43, status='unknown', file='output/out_hourly_tower.csv')
          write(43,9190) 'iyear', 'jday', 'istep-1', 'swin', 'swout', 'pari', 'paro',  &
                                  'apar', 'rn','-fsena(1)','-fvapa(1)*hvap','-soihfl(1)',   & 
                                  '-tneetot(1)*1e6','plai(1,16)'
9190  format (13(a,x))

          open (41, status='unknown', file='output/out_daily_tower.csv')
	  write(41,9192)'iyear','jday','flx(1)','flx(2)','flx(3)','flx(4)',          &
                                'flx(5)','flx(6)','flx(7)','flx(8)','flx(9)','flx(11)','flx(10)',  &
                                'plai(1,16)*grnfraccrop(1,16)','plai(1,16)','flx(12)',     &
                                'flx(13)','flx(14)','flx(15)'
9192  format (18(a,x))

         open (225, status='unknown', file='output/soil.csv')
	 write(225,9193)'jday','wsoi(1,1)','wsoi(1,2)','wsoi(1,3)','wsoi(1,4)','wsoi(1,5)','wsoi(1,6)',  &
                         'wsoi(1,7)','wsoi(1,8)','wsoi(1,9)','wsoi(1,10)','wsoi(1,11)',               &
                         'sfield(1,1)','swilt(1,1)','stresstl(1)'
9193  format (15(a,x))

         open (42, status='unknown', file='output/crops_out.csv')
	 write(42,9194) 'iyear','jday', 'idpp(i,j)','hui(i,j)', 'aroot(i,j)',&
                                     'aleaf(i,j)', 'astem(j)','arepr(j)','plai(i,j)'  
9194  format (9(a,x))

         open (222, status='unknown', file='output/biomass.csv')
	 write(222,9195) 'iyear','jday','idpp','plai','plai*grnfraccrop',             &
                         'cbiol*(1/0.45)*10.0','(aylprod-cbiol)*(1/0.45)*10.0',      &
                         '(aybprod-ayrprod-aylprod-cbiog)*(1/0.45)*10.0',     &
                         'cbiog*(1/0.45)*10.0,(aybprod-ayrprod)*(1/0.45)*10.0',    &
                         'ayrprod*(1/0.45)*10.0' 
9195  format (10(a,x))

         open (23, status='unknown', file='output/physiology.csv')
	 write(23,9196) 'jday','dtime/3600.','lai(i,1)*fl(i)*grnfraccrop',               &
                        'js*1e+06','je*1e+06','jc*1e+06','anc3*1e+06','tleaf','tempvm'
9196  format (9(a,x))

         open (226, status='unknown', file='output/mdaily.csv')
	 write(226,9197) 'iyear','jday', 'adrain'

9197  format (3(a,x))

         open (230, status='unknown', file='output/water.csv')
	 write(230,9198) 'iyear', 'jday', 'adrain' ,'adaet','adevap','adtrans'          &
                         ,'adwsoi*poros(1,1)*20.*10.', 'adwsoi2*poros(1,1)*30.*10.'    &
                         ,'(adwsoi*0.4 + adwsoi2*0.6)*poros(1,1)*50.*10.'            &
                         ,'adtrunoff','adsrunoff','addrainage'

9198  format (12(a,x))

#endif /* SINGLE_POINT_MODEL */

      return
end subroutine build_file
