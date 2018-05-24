# include <stdlib.h>

# include "../optis_global.h"

paramsrange_opt* alloc_paramsrange(int ncal, calibration_lvl **allcaliblvl)
{
	paramsrange_opt *prange;
	int i;

	prange = (paramsrange_opt*)malloc(sizeof(paramsrange_opt));

	prange->ncal = ncal;

	prange->min_value = (double**)malloc(ncal*sizeof(double*));
	prange->max_value = (double**)malloc(ncal*sizeof(double*));

	for(i=0;i<ncal;i++)
	{
		prange->min_value[i] = (double*)malloc(allcaliblvl[i]->cpars*sizeof(double));
		prange->max_value[i] = (double*)malloc(allcaliblvl[i]->cpars*sizeof(double));
	}

	return prange;
}
