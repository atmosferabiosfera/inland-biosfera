# include <stdlib.h>
# include <unistd.h>

# include "optis_global.h"

int main()
{
    int ncal, i;
    population *pop;
    individual *best_ind;
    nsga_param *nparam;
    modelo_st* modelo;
    calibration_lvl *calib, *allobj;
    calibration_lvl **allcaliblvl;
    base_param *bparam;
    sit_array *sitios;
    FILE *f_bestind, *f_calibparams, *f_conf;
    char conf_name[100], cmd[200];
    paramsrange_opt *prange_opt;

    /* Limpando o diretorio de output */
    sprintf(cmd, "rm -rf %s*;", OUTPUTDIR);
    system(cmd);

    sprintf(conf_name, "%s%s", INPUTDIR, CONF_FILE);
    if(access(conf_name, R_OK) != 0){
    	printf("\nERRO: O arquivo %s nao existe ou nao pode ser lido.\n", conf_name);
    	exit(EXIT_FAILURE);
    }
    f_conf = fopen(conf_name, "r");

    /* Ler os parametros iniciais */
    nparam = init_nsgaparam(f_conf);
    modelo = init_modelo(f_conf);
    bparam = init_baseparam(f_conf);
    sitios = init_sitarray(f_conf);

    create_tempdir(nparam->popsize, sitios, modelo);

    f_bestind = open_bestindf();
    f_calibparams = open_calibparamsf();

    allobj = alloc_allobj(sitios);

    /* Ler o total de calibracoes (ncal) */
    ncal = rd_int(f_conf);
    if(ncal < 1){
    	printf("\nERRO: O total de niveis hierarquicos precisa ser maior que zero.\n");
    	exit(EXIT_FAILURE);
    }

    printf("%d calibration levels will be performed.\n", ncal);

    allcaliblvl = init_allcalibrationlvl(f_conf, bparam, sitios, ncal, allobj);

    prange_opt = alloc_paramsrange(ncal, allcaliblvl);

    for(i=0;i<ncal;i++)
    {
    	printf("Level %d\n", i+1);

    	/* Ler os dados relativos a calibracao */
    	calib = allcaliblvl[i];

    	pop = nsga(nparam, calib, bparam, sitios, modelo);

    	report_paretofront(bparam, nparam->popsize, pop, calib);

    	set_paramsrange(i, nparam->popsize, pop, calib->cpars, prange_opt);

    	best_ind = mindist_origin(pop, nparam->popsize, calib->nobj);

    	store_params(best_ind, calib, nparam->popsize, bparam, sitios);

    	report_bestind(f_bestind, best_ind, calib);
    	report_calibparams(calib, f_calibparams, best_ind, calib->parname);

    	dealloc_pop(pop, nparam->popsize);
    }

    report_inlandparams(bparam);

    report_allobj(allobj, bparam, sitios, modelo);

    report_paramsrange(prange_opt, allcaliblvl);

    fclose(f_bestind);
    fclose(f_calibparams);
    fclose(f_conf);

    dealloc_allcalibrationlvl(ncal, allcaliblvl, sitios);
    dealloc_calibrationlvl(allobj, sitios);
    dealloc_paramsrange(prange_opt);

    dealloc_nsgaparam(nparam);
    dealloc_modelo(modelo);
    dealloc_baseparam(bparam);
    dealloc_sitarray(sitios);

    printf("End of calibration\n");

    return EXIT_SUCCESS;
}
