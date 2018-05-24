# include <stdlib.h>

# include "optis_global.h"

individual* mindist_origin(population *pop, int popsize, int nobj)
{
	double *objmax, dmin, d;
	individual *best_ind;
	int i, j;

	objmax = (double*)malloc(nobj*sizeof(double));

	for(i=0;i<nobj;i++)
	{
		objmax[i] = pop->ind[0].obj[i];
	}
	for(i=1;i<popsize;i++)
	{
		for(j=0;j<nobj;j++)
		{
			if(pop->ind[i].obj[j] > objmax[j])
			{
				objmax[j] = pop->ind[i].obj[j];
			}
		}
	}

	dmin = calcdist_origin(&(pop->ind[0]), objmax, nobj);
	best_ind = &(pop->ind[0]);

	for(i=1;i<popsize;i++)
	{
		d = calcdist_origin(&(pop->ind[i]), objmax, nobj);
		if(d < dmin)
		{
			dmin = d;
			best_ind = &(pop->ind[i]);
		}
	}

	free(objmax);

	return best_ind;
}
