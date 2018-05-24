# include <stdlib.h>

# include "optis_global.h"

void store_params(individual *ind, calibration_lvl *calib, int popsize, base_param *bparams, sit_array *sitios)
{
    int i, j;
    char *path;
    constraint_array *constraint;

    constraint = calib->constraint;

    path = (char*)malloc(500*sizeof(char));

    for(i=0;i<popsize;i++)
    {
    	for(j=0;j<sitios->size;j++)
    	{
    		sprintf(path, "%s%d/%s%d/", TEMPDIR, i, B_SITIODIR, j);

    		write_inlandparams(path, ind->xreal, bparams, calib);
    	}
    }

    // Diretório auxiliar
    sprintf(path, "%s%s/", TEMPDIR, AUXBASEDIR);
    write_inlandparams(path, ind->xreal, bparams, calib);


    free(path);

//    if(constraint->size > 0)
//    {
//    	for (i = 0; i < constraint->size; ++i) {
//    		stparams[constraint->cons[i][3]] = 1.0 - (stparams[constraint->cons[i][1]] + stparams[constraint->cons[i][2]]);
//    	}
//    }

    return;
}
