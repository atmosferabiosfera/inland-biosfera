# include <omp.h>

# include "optis_global.h"

void evaluate_pop(population *pop, int popsize, calibration_lvl *calib, base_param *bparam, sit_array *sitios, modelo_st *modelo)
{
    int i;

    #pragma omp parallel private(i)
    {
	    #pragma omp for
    	for(i=0;i<popsize;i++)
    	{
    		printf("Evaluating ind %d\n", i);
    		evaluate_ind(i, &(pop->ind[i]), calib, bparam, sitios, modelo);
    	}
    }

    return;
}
