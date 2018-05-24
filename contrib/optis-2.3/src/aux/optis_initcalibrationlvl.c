# include <stdlib.h>
# include <string.h>

# include "../optis_global.h"


calibration_lvl* init_calibrationlvl(FILE *file, int lvl, base_param *bparam, sit_array *sitios)
{
    calibration_lvl *calib;
    int i, j;

    calib = (calibration_lvl*)malloc(sizeof(calibration_lvl));

    calib->obj = (var_array**)malloc(sitios->size*sizeof(var_array*));
    for (j = 0; j < sitios->size; j++) {
		calib->obj[j] = NULL;
	}

    calib->constraint = alloc_constraint();

    calib->lvl = lvl;

    //Lendo o modo de calibracao
    calib->mode = rd_int(file);

    if((calib->mode != CALIB_MOD_IND) && (calib->mode != CALIB_MOD_MED)){
    	printf("\nERRO: Esse modo de calibracao nao existe.\n");
    	exit(EXIT_FAILURE);
    }

    calib->name = rd_string(file, 50);

    /* Lendo os parametros */
    calib->cpars = rd_int(file);
    if(calib->cpars < 1){
    	printf("\nERRO: O total de parametros calibrados precisa ser maior que zero.\n");
    	exit(EXIT_FAILURE);
    }

    calib->parl = (int*)malloc(calib->cpars*sizeof(int));
    calib->parc = (int*)malloc(calib->cpars*sizeof(int));
    calib->pars = (int*)malloc(calib->cpars*sizeof(int));
    calib->parf = (int*)malloc(calib->cpars*sizeof(int));
    calib->parname = (char**)malloc(calib->cpars*sizeof(char*));
    calib->min_parvalue = (double*)malloc(calib->cpars*sizeof(double));
    calib->max_parvalue = (double*)malloc(calib->cpars*sizeof(double));

    for(i=0;i<calib->cpars;i++)
    {
    	calib->parname[i] = rd_string(file, 50);

    	calib->parf[i] = rd_int(file)-1;
    	calib->parl[i] = rd_int(file);
    	calib->parc[i] = rd_int(file);
    	calib->pars[i] = rd_int(file);

    	if(calib->pars[i] < 7){
    		printf("\nERRO no parametro %s: O tamanho do parametro precisa ser maior que 6. Tamanho definido: %d\n", calib->parname[i], calib->pars[i]);
    		exit(EXIT_FAILURE);
    	}


    	calib->min_parvalue[i] = rd_double(file);
    	calib->max_parvalue[i] = rd_double(file);

    	 if(calib->min_parvalue[i] >= calib->max_parvalue[i]){
    		 printf("\nERRO: No parametro %d o limite minimo precisa ser menor que o maximo.\n", i+1);
    		 exit(EXIT_FAILURE);
    	 }

//    	check_paramconstraint(param_name, calib->parpos[i], calib->constraint);

    }

//    check_constraints(calib->constraint, bparam);

    /* Lendo as funcoes objetivas */
	rd_calibrationvars(lvl, calib, file, sitios);

    return calib;
}
