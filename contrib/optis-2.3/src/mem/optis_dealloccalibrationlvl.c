# include <stdlib.h>

# include "../optis_global.h"

void dealloc_calibrationlvl(calibration_lvl *calib, sit_array *sitios)
{
    int i, k;

    free(calib->parf);
    free(calib->parl);
    free(calib->parc);
    free(calib->pars);
    free(calib->min_parvalue);
    free(calib->max_parvalue);
    free(calib->name);
    free(calib->obj_comp);

    for(i=0;i<calib->cpars;i++)
    {
    	free(calib->parname[i]);
    }
    free(calib->parname);

    for(i=0;i<calib->nobj;i++)
    {
    	free(calib->of_name[i]);
    }
    free(calib->of_name);

    for(i=0;i<sitios->size;i++)
    {
    	if(calib->obj[i] != NULL)
    	{
    		for(k=0;k<N_OF_METHODS;k++)
    		{
    			if(calib->obj[i][k].size != 0)
    			{
    				free(calib->obj[i][k].var);
    			}
    		}
    		free(calib->obj[i]);
    	}
    }
    free(calib->obj);

    dealloc_constraint(calib->constraint);

    free(calib);

    return;
}
