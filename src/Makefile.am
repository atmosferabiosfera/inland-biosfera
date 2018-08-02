## Process this file with automake to produce Makefile.in
## In turn, Makefile.in will be used by the configure script to generate the
## final Makefile.

if SINGLE_POINT

# We compile a different set of files if making a single point model

bin_PROGRAMS = inland-single_point

 inland_single_point_SOURCES = inland_main_offline.F90 inland_arr2vec.F90 inland_canini.F90 inland_canopy.F90 inland_cascad2.F90 inland_cascade.F90 inland_co2.F90 inland_coldstart.F90 inland_com1d.F90 inland_comatm.F90 inland_combcs.F90 inland_combgc.F90 inland_comdiag.F90 inland_comfire.F90 inland_comforc.F90 inland_comcrop.F90 inland_comhyd.F90 inland_comhour.F90 inland_commenthandler.F90 inland_compft.F90 inland_comsatp.F90 inland_comsno.F90 inland_comsoi.F90 inland_comsum.F90 inland_comtex.F90 inland_comveg.F90 inland_comwork.F90 inland_const.F90 inland_daily.F90 inland_delta_sw.F90 inland_diurnal.F90 inland_drystress.F90 inland_dynaveg.F90 inland_existence.F90 inland_firectem.F90 inland_fstrat.F90 inland_fwetcal.F90 inland_gdiag.F90 inland_alloc.F90 inland_control.F90 inland_map_offline.F90 inland_parameters.F90 inland_initsw.F90 inland_impexp.F90 inland_impexp2.F90 inland_inisnow.F90 inland_inisoil.F90 inland_inisum.F90 inland_inisurf.F90 inland_initialib.F90 inland_iniveg.F90 inland_inpoly.F90 inland_irrad.F90 inland_lenchrib.F90 inland_linsolve.F90 inland_lsmmapib.F90 inland_lsxmain.F90 inland_mix.F90 inland_noveg.F90 inland_pheno.F90 inland_prealloc.F90 inland_inneralloc.F90 inland_random.F90 inland_rdparam.F90 inland_rdenv.F90 inland_readforc.F90 inland_readit.F90 inland_readit_sparse.F90 inland_readitem.F90 inland_readitems.F90 inland_scopya.F90 inland_setarrib.F90 inland_setsoi.F90 inland_snow.F90 inland_snowheat.F90 inland_soilbgc.F90 inland_soilctl.F90 inland_soilh2o.F90 inland_soilheat.F90 inland_solalb.F90 inland_solarf.F90 inland_solset.F90 inland_solsur.F90 inland_steph2o.F90 inland_steph2o2.F90 inland_stomataib.F90 inland_subgrid.F90 inland_sumday.F90 inland_summonth.F90 inland_sumnow.F90 inland_sumyear.F90 inland_textcls.F90 inland_tridiaib.F90 inland_tscreen.F90 inland_turcof.F90 inland_turvap.F90 inland_twet3.F90 inland_twoset.F90 inland_twostrib.F90 inland_vadapt.F90 inland_vdiag.F90 inland_vec2arr.F90 inland_vegmap.F90 inland_wadjust.F90 inland_wdiag.F90 inland_comnitr.F90 inland_nitrostress.F90 inland_phenocrop.F90 inland_planting.F90 inland_rotation.F90 inland_vernalization.F90 inland_cropback.F90 inland_cropoutput.F90 inland_cropresidue.F90 inland_cropupdate.F90 inland_initialcrop.F90 inland_irrigation.F90 inland_leaching.F90 inland_constint.F90 inland_dailymet.F90 inland_single.F90 inland_build_file.F90 inland_check.F90

else

# For the 2D ordinary model, we use these files

bin_PROGRAMS = inland-grid

 inland_grid_SOURCES = inland_main_offline.F90 inland_arr2vec.F90 inland_canini.F90 inland_canopy.F90 inland_cascad2.F90 inland_cascade.F90 inland_climanl.F90 inland_climanl2.F90 inland_co2.F90 inland_coldstart.F90 inland_com1d.F90 inland_comatm.F90 inland_combcs.F90 inland_combgc.F90 inland_comcrop.F90 inland_comdiag.F90 inland_comfire.F90 inland_comhour.F90  inland_comhyd.F90 inland_commenthandler.F90 inland_compft.F90 inland_comsatp.F90 inland_comsno.F90 inland_comsoi.F90 inland_comsum.F90 inland_comtex.F90 inland_comveg.F90 inland_comwork.F90 inland_const.F90 inland_daily.F90 inland_dailymon.F90 inland_delta_sw.F90 inland_diurnal.F90 inland_drystress.F90 inland_dynaveg.F90 inland_existence.F90 inland_firectem.F90 inland_fstrat.F90 inland_fwetcal.F90 inland_gdiag.F90 inland_alloc.F90 inland_control.F90 inland_map_offline.F90 inland_parameters.F90 inland_initsw.F90 inland_impexp.F90 inland_impexp2.F90 inland_inifile.F90 inland_inird.F90 inland_inisnow.F90 inland_inisoil.F90 inland_inisum.F90 inland_inisurf.F90 inland_initialib.F90 inland_inivar.F90 inland_iniveg.F90 inland_inpoly.F90 inland_irrad.F90 inland_lenchrib.F90 inland_linsolve.F90 inland_lsmmapib.F90 inland_lsxmain.F90 inland_mix.F90 inland_netcdfutils.F90 inland_noveg.F90 inland_pheno.F90 inland_prealloc.F90 inland_inneralloc.F90 inland_random.F90 inland_rdday.F90 inland_rdlndtrans.F90 inland_rdmon.F90 inland_rdparam.F90 inland_rdenv.F90 inland_readit.F90 inland_readit_sparse.F90 inland_readitem.F90 inland_readitems.F90 inland_readvar.F90 inland_restartib.F90 inland_scopya.F90 inland_setarrib.F90 inland_setsoi.F90 inland_snow.F90 inland_snowheat.F90 inland_soilbgc.F90 inland_soilctl.F90 inland_soilh2o.F90 inland_soilheat.F90 inland_solalb.F90 inland_solarf.F90 inland_solset.F90 inland_solsur.F90 inland_steph2o.F90 inland_steph2o2.F90 inland_stomataib.F90 inland_subgrid.F90 inland_sumday.F90 inland_summonth.F90 inland_sumnow.F90 inland_sumyear.F90 inland_textcls.F90 inland_tridiaib.F90 inland_tscreen.F90 inland_turcof.F90 inland_turvap.F90 inland_twet3.F90 inland_twoset.F90 inland_twostrib.F90 inland_vadapt.F90 inland_vdiag.F90 inland_vec2arr.F90 inland_vegmap.F90 inland_wadjust.F90 inland_wdaily.F90 inland_wdiag.F90 inland_wmonthly.F90 inland_wrestart.F90 inland_writevar.F90 inland_wyearly.F90 inland_comnitr.F90 inland_nitrostress.F90 inland_phenocrop.F90 inland_planting.F90 inland_rotation.F90 inland_vernalization.F90 inland_cropback.F90 inland_cropoutput.F90 inland_cropresidue.F90 inland_cropupdate.F90 inland_initialcrop.F90 inland_irrigation.F90 inland_leaching.F90 inland_constint.F90 inland_diurnalmet.F90 inland_dailymet.F90 inland_methourly.F90 inland_arr3arr2.F90 inland_logicf.F90 inland_logict.F90 inland_rdhrvegmap.F90 inland_whrvegmap.F90 inland_single_agro.F90 inland_build_file.F90 inland_check.F90

endif

CLEANFILES=*.mod *.oo

inland_arr2vec.o: inland_arr2vec.F90 inland_parameters.o inland_combcs.o inland_subgrid.o
inland_canopy.o: inland_canopy.F90 inland_com1d.o inland_parameters.o inland_comveg.o inland_comatm.o 
inland_cascad2.o: inland_cascad2.F90 inland_com1d.o inland_parameters.o inland_comveg.o inland_comatm.o inland_comsno.o 
inland_cascade.o: inland_cascade.F90 inland_com1d.o inland_parameters.o inland_comveg.o inland_comatm.o 
inland_climanl.o: inland_climanl.F90 inland_parameters.o inland_comveg.o inland_combcs.o inland_comcrop.o
inland_climanl2.o: inland_climanl2.F90 inland_parameters.o inland_comveg.o inland_combcs.o 
inland_co2.o: inland_co2.F90 inland_control.o
inland_coldstart.o: inland_coldstart.F90 inland_parameters.o inland_comsoi.o inland_comsno.o inland_combcs.o 
inland_com1d.o: inland_com1d.F90 
inland_comatm.o: inland_comatm.F90 
inland_combcs.o: inland_combcs.F90 
inland_combgc.o: inland_combgc.F90
inland_comcrop.o: inland_comcrop.F90 
inland_comdiag.o: inland_comdiag.F90 
inland_comfire.o: inland_comfire.F90
inland_comhyd.o: inland_comhyd.F90 
inland_comhour.o: inland_comhour.F90
inland_commenthandler.o: inland_commenthandler.F90 
inland_compft.o: inland_compft.F90 
inland_comsatp.o: inland_comsatp.F90 
inland_comsno.o: inland_comsno.F90 
inland_comsoi.o: inland_comsoi.F90 
inland_comsum.o: inland_comsum.F90 
inland_comtex.o: inland_comtex.F90 
inland_comveg.o: inland_comveg.F90 
inland_comwork.o: inland_comwork.F90 inland_parameters.o 
inland_comnitr.o: inland_comnitr.F90
inland_const.o: inland_const.F90 
inland_check.o: inland_check.F90 inland_control.o inland_comcrop.F90
inland_daily.o: inland_daily.F90 ../include/inland_comsat.h inland_parameters.o inland_combcs.o inland_comatm.o inland_comsum.o inland_comveg.o inland_comsatp.o inland_control.o inland_subgrid.o
inland_diurnal.o: inland_diurnal.F90 ../include/inland_comsat.h inland_parameters.o inland_comatm.o inland_comwork.o inland_comsatp.o inland_comnitr.o  inland_subgrid.o
inland_drystress.o: inland_drystress.F90 inland_parameters.o inland_comveg.o inland_comsoi.o 
inland_existence.o: inland_existence.F90 inland_parameters.o inland_compft.o inland_comveg.o inland_comnitr.o
inland_firectem.o: inland_firectem.F90 inland_comveg.o inland_comsoi.o inland_comfire.o inland_comatm.o inland_comwork.o inland_random.o
inland_fstrat.o: inland_fstrat.F90 inland_parameters.o 
inland_fwetcal.o: inland_fwetcal.F90 inland_com1d.o inland_parameters.o inland_comveg.o 
inland_gdiag.o: inland_gdiag.F90 inland_parameters.o inland_comveg.o inland_comsum.o inland_comnitr.o inland_control.o inland_subgrid.o
inland_control.o: inland_control.F90 
inland_map_offline.o: inland_map_offline.F90 inland_parameters.o inland_comwork.o inland_lsmmapib.o 
inland_parameters.o: inland_parameters.F90
inland_impexp.o: inland_impexp.F90 inland_parameters.o 
inland_impexp2.o: inland_impexp2.F90 inland_parameters.o 
inland_inird.o: inland_inird.F90 
inland_inisnow.o: inland_inisnow.F90 inland_parameters.o inland_comsoi.o inland_comsno.o 
inland_inisoil.o: inland_inisoil.F90 inland_parameters.o inland_comsoi.o inland_comsno.o inland_comtex.o inland_comsum.o inland_comcrop.o inland_comnitr.o
inland_inisum.o: inland_inisum.F90 inland_parameters.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o 
inland_inisurf.o: inland_inisurf.F90 inland_parameters.o inland_comveg.o inland_comsoi.o inland_comatm.o inland_comhyd.o 
inland_initialib.o: inland_initialib.F90 inland_parameters.o 
inland_iniveg.o: inland_iniveg.F90 inland_parameters.o inland_combcs.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_compft.o 
inland_inpoly.o: inland_inpoly.F90 
inland_irrad.o: inland_irrad.F90 inland_parameters.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_lenchrib.o: inland_lenchrib.F90 
inland_linsolve.o: inland_linsolve.F90 inland_parameters.o
inland_lsmmapib.o: inland_lsmmapib.F90 
inland_lsxmain.o: inland_lsxmain.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_com1d.o 
inland_mix.o: inland_mix.F90 inland_parameters.o 
inland_noveg.o: inland_noveg.F90 inland_parameters.o inland_comsno.o inland_comsoi.o inland_comveg.o 
inland_pheno.o: inland_pheno.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comcrop.o
inland_prealloc.o: inland_prealloc.F90 inland_parameters.o inland_compft.o inland_comveg.o inland_comtex.o inland_combgc.o inland_comsoi.o inland_comwork.o
inland_inneralloc.o: inland_inneralloc.F90 inland_parameters.o inland_comsoi.o inland_comwork.o
inland_random.o: inland_random.F90 inland_control.o
inland_rdparam.o: inland_rdparam.F90 inland_parameters.o inland_comveg.o inland_comsoi.o inland_compft.o inland_comtex.o inland_combgc.o 
inland_rdenv.o: inland_rdenv.F90 inland_parameters.o inland_control.o
inland_readitem.o: inland_readitem.F90 
inland_readitems.o: inland_readitems.F90 
inland_scopya.o: inland_scopya.F90 
inland_setarrib.o: inland_setarrib.F90 
inland_setsoi.o: inland_setsoi.F90 ../include/inland_comsat.h inland_parameters.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_comsatp.o 
inland_snow.o: inland_snow.F90 inland_parameters.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_snowheat.o: inland_snowheat.F90 inland_parameters.o inland_comsno.o inland_comsoi.o 
inland_soilbgc.o: inland_soilbgc.F90 inland_parameters.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_control.o
inland_soilctl.o: inland_soilctl.F90 inland_parameters.o inland_control.o inland_comatm.o inland_comhyd.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_soilh2o.o: inland_soilh2o.F90 inland_parameters.o inland_comsoi.o inland_com1d.o 
inland_soilheat.o: inland_soilheat.F90 inland_parameters.o inland_comsoi.o 
inland_solalb.o: inland_solalb.F90 inland_parameters.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_solarf.o: inland_solarf.F90 inland_parameters.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_solset.o: inland_solset.F90 inland_parameters.o inland_comatm.o inland_comsno.o inland_comveg.o inland_com1d.o 
inland_solsur.o: inland_solsur.F90 inland_parameters.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_steph2o.o: inland_steph2o.F90 inland_parameters.o 
inland_steph2o2.o: inland_steph2o2.F90 ../include/inland_comsat.h inland_parameters.o inland_comsatp.o 
inland_stomataib.o: inland_stomataib.F90 ../include/inland_comsat.h inland_parameters.o inland_comatm.o inland_comsatp.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_com1d.o inland_compft.o 
inland_subgrid.o: inland_subgrid.F90 inland_parameters.o inland_comwork.o 
inland_sumday.o: inland_sumday.F90 inland_parameters.o inland_comatm.o inland_comhyd.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_com1d.o 
inland_summonth.o: inland_summonth.F90 inland_parameters.o inland_comatm.o inland_comhyd.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_com1d.o inland_control.o
inland_sumnow.o: inland_sumnow.F90 inland_parameters.o inland_control.o inland_comatm.o inland_comhyd.o inland_comsoi.o inland_comsno.o inland_comsum.o inland_comveg.o inland_com1d.o 
inland_sumyear.o: inland_sumyear.F90 inland_com1d.o inland_parameters.o inland_comatm.o inland_comhyd.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_control.o
inland_textcls.o: inland_textcls.F90
inland_tridiaib.o: inland_tridiaib.F90 inland_parameters.o 
inland_tscreen.o: inland_tscreen.F90 inland_parameters.o 
inland_turcof.o: inland_turcof.F90 inland_parameters.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_turvap.o: inland_turvap.F90 ../include/inland_comsat.h inland_parameters.o inland_control.o inland_comatm.o inland_comhyd.o inland_comsatp.o inland_comsno.o inland_comsoi.o inland_comveg.o inland_com1d.o 
inland_twet3.o: inland_twet3.F90 ../include/inland_comsat.h inland_parameters.o inland_comsatp.o 
inland_twoset.o: inland_twoset.F90 inland_parameters.o inland_comveg.o inland_com1d.o 
inland_twostrib.o: inland_twostrib.F90 inland_parameters.o inland_comveg.o inland_com1d.o 
inland_vadapt.o: inland_vadapt.F90 inland_parameters.o 
inland_vdiag.o: inland_vdiag.F90 inland_parameters.o inland_comveg.o inland_comsum.o inland_control.o inland_subgrid.o
inland_vec2arr.o: inland_vec2arr.F90 inland_parameters.o inland_comwork.o inland_subgrid.o
inland_vegmap.o: inland_vegmap.F90 inland_parameters.o inland_comveg.o inland_subgrid.o
inland_wadjust.o: inland_wadjust.F90 inland_comsoi.o inland_com1d.o inland_parameters.o inland_comhyd.o 
inland_wdiag.o: inland_wdiag.F90 inland_parameters.o inland_com1d.o inland_comveg.o inland_comhyd.o inland_comatm.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comdiag.o inland_comnitr.o inland_control.o
inland_rdhrvegmap.o: inland_comwork.o inland_subgrid.o inland_combcs.o
inland_delta_sw.o: inland_delta_sw.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comveg.o
inland_initsw.o:  inland_delta_sw.F90 inland_parameters.o inland_comatm.o inland_comsoi.o 


if SINGLE_POINT

inland_canini.o: inland_canini.F90 inland_com1d.o inland_parameters.o inland_control.o inland_comveg.o inland_comsoi.o inland_comatm.o inland_comsno.o 
inland_comforc.o: inland_comforc.F90
inland_alloc.o: inland_alloc.F90 inland_com1d.o inland_lsmmapib.o inland_parameters.o inland_control.o inland_comatm.o inland_combcs.o inland_comhyd.o inland_comfire.o inland_compft.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comforc.o inland_comwork.o
inland_dynaveg.o: inland_dynaveg.F90 inland_parameters.o inland_compft.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_control.o inland_comfire.o
inland_main_offline.o: inland_main_offline.F90 inland_parameters.o inland_comatm.o inland_combgc.o inland_comcrop.o inland_comdiag.o inland_comwork.o inland_lsmmapib.o inland_control.o inland_comhour.o inland_comveg.o inland_com1d.o inland_comsoi.o inland_comhyd.o inland_comsum.o inland_comforc.o inland_rdenv.o inland_comfire.o inland_single.o
inland_readforc.o: inland_readforc.F90 inland_parameters.o inland_control.o inland_combcs.o inland_comsoi.o inland_comveg.o inland_comforc.o inland_parameters.o
inland_readit.o: inland_readit.F90 inland_parameters.o inland_combcs.o inland_comsoi.o inland_comveg.o inland_comwork.o inland_subgrid.o inland_control.o
inland_readit_sparse.o: inland_readit_sparse.F90 inland_parameters.o inland_combcs.o inland_comsoi.o inland_comveg.o inland_comwork.o inland_subgrid.o inland_control.o

inland_single.o: inland_single.F90 inland_comatm.o inland_comsoi.o inland_com1d.o inland_comveg.o inland_control.o inland_parameters.o inland_comhyd.o inland_comsum.o inland_comforc.o
inland_build_file.o: inland_build_file.F90 inland_comatm.o inland_com1d.o inland_comcrop.o inland_comveg.o inland_control.o inland_comsoi.o inland_comhyd.o

else # if 2D Model

inland_canini.o: inland_canini.F90 inland_com1d.o inland_parameters.o inland_comveg.o inland_comsoi.o inland_comatm.o inland_comsno.o 
inland_dailymon.o: inland_dailymon.F90 ../include/inland_comsat.h inland_parameters.o inland_combcs.o inland_comatm.o inland_comsum.o inland_comveg.o inland_comsatp.o inland_control.o inland_combcs.o inland_subgrid.o
inland_alloc.o: inland_alloc.F90 inland_com1d.o inland_lsmmapib.o inland_parameters.o inland_control.o inland_comatm.o inland_combcs.o inland_comhyd.o inland_comfire.o inland_compft.o inland_comsno.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comwork.o
inland_dynaveg.o: inland_dynaveg.F90 inland_parameters.o inland_compft.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_control.o inland_comfire.o inland_rdlndtrans.o
inland_inifile.o: inland_inifile.F90 inland_control.o
inland_inivar.o: inland_inivar.F90 inland_comwork.o inland_control.o
inland_main_offline.o: inland_main_offline.F90 inland_parameters.o inland_comatm.o inland_combgc.o inland_comcrop.o inland_comdiag.o inland_comwork.o inland_lsmmapib.o inland_control.o inland_rdenv.o inland_comfire.o inland_comhour.o inland_comveg.o inland_com1d.o inland_comhyd.o inland_comsoi.o inland_random.o inland_subgrid.o inland_whrvegmap.o
inland_netcdfutils.o: inland_netcdfutils.F90
inland_rdday.o: inland_rdday.F90 inland_parameters.o inland_combcs.o inland_comwork.o inland_control.o
inland_rdlndtrans.o: inland_rdlndtrans.F90 inland_parameters.o inland_combcs.o inland_control.o inland_comveg.o inland_comwork.o
inland_rdmon.o: inland_rdmon.F90 inland_parameters.o inland_combcs.o inland_comwork.o inland_control.o
inland_readit.o: inland_readit.F90 inland_parameters.o inland_combcs.o inland_comsoi.o inland_comveg.o inland_comwork.o inland_subgrid.o inland_control.o 
inland_readit_sparse.o: inland_readit_sparse.F90 inland_parameters.o inland_combcs.o inland_comsoi.o inland_comveg.o inland_comwork.o inland_subgrid.o inland_control.o 
inland_rdhrvegmap.o
inland_readvar.o: inland_readvar.F90 inland_comwork.o
inland_restartib.o: inland_restartib.F90 inland_parameters.o inland_control.o inland_comveg.o inland_comsoi.o inland_comsno.o inland_combcs.o inland_comwork.o inland_comsum.o inland_comatm.o
inland_wdaily.o: inland_wdaily.F90 inland_parameters.o inland_comatm.o inland_comwork.o inland_comsum.o inland_control.o inland_netcdfutils.o inland_subgrid.o
inland_wmonthly.o: inland_wmonthly.F90 inland_parameters.o inland_comsoi.o inland_comwork.o inland_comsum.o inland_control.o inland_netcdfutils.o inland_subgrid.o
inland_wrestart.o: inland_wrestart.F90 inland_parameters.o inland_comveg.o inland_comsoi.o inland_comsno.o inland_comwork.o inland_comsum.o inland_control.o inland_comatm.o inland_combcs.o inland_netcdfutils.o
inland_writevar.o: inland_writevar.F90 inland_comwork.o inland_parameters.o inland_subgrid.o inland_netcdfutils.o
inland_wyearly.o: inland_wyearly.F90 inland_parameters.o inland_comveg.o inland_comwork.o inland_comsum.o inland_control.o inland_subgrid.o inland_netcdfutils.o inland_comfire.o
inland_whrvegmap.o: inland_whrvegmap.F90 inland_parameters.o inland_control.o inland_combcs.o inland_subgrid.o inland_netcdfutils.o
inland_cropoutput.o: inland_cropoutput.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_cropresidue.o: inland_cropresidue.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_cropupdate.o: inland_cropupdate.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsno.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o inland_compft.o
inland_initialcrop.o: inland_initialcrop.F90 inland_parameters.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_irrigation.o: inland_irrigation.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_leaching.o: inland_leaching.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_nitrostress.o: inland_nitrostress.F90 inland_parameters.o inland_com1d.o inland_comatm.o inland_comsoi.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_phenocrop.o: inland_phenocrop.F90  inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsno.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o inland_compft.o
inland_planting.o: inland_planting.F90 inland_phenocrop.F90  inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsno.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o inland_combcs.o
inland_rotation.o: inland_rotation.F90 inland_parameters.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_vernalization.o: inland_vernalization.F90 inland_parameters.o inland_comatm.o inland_comsoi.o inland_comsno.o inland_comsum.o inland_comveg.o inland_comcrop.o inland_comnitr.o
inland_cropback.o: inland_cropback.F90 inland_parameters.o
inland_constint.o: inland_constint.F90
inland_diurnalmet.o: inland_diurnalmet.F90 inland_comveg.o inland_comhour.o inland_comcrop.o inland_comnitr.o inland_parameters.o inland_comatm.o inland_comwork.o inland_comsatp.o ../include/inland_comsat.h
inland_dailymet.o: inland_dailymet.F90 inland_parameters.o inland_control.o inland_comveg.o inland_comhour.o inland_combcs.o inland_comatm.o inland_comsoi.o inland_comcrop.o inland_comsatp.o ../include/inland_comsat.h
inland_methourly.o: inland_methourly.F90 inland_comwork.o inland_comhyd.o inland_comsum.o inland_comhour.o inland_comnitr.o ../include/inland_comsat.h
inland_arr3arr2.o: inland_arr3arr2.F90 inland_parameters.o
inland_logicf.o: inland_logicf.F90
inland_single_agro.o: inland_single_agro.F90 inland_comatm.o inland_com1d.o inland_comcrop.o inland_comveg.o inland_control.o inland_comsoi.o inland_comhyd.o
inland_build_file.o: inland_build_file.F90 inland_comatm.o inland_com1d.o inland_comcrop.o inland_comveg.o inland_control.o inland_comsoi.o inland_comhyd.o
endif


#	@echo Compiling $< into $@ \
#	$(F90) -c $(F90FLAGS) $(CPPFLAGS) $< -o $@
