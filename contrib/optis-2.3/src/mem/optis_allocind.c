# include <stdlib.h>

# include "../optis_global.h"

void alloc_ind(individual *ind, int cpars, int nobj)
{
	if(cpars > 0)
	{
		ind->xreal = (double*)malloc(cpars*sizeof(double));
	}
	else
	{
		ind->xreal = NULL;
	}

	if(nobj > 0)
	{
		ind->obj = (double*)malloc(nobj*sizeof(double));
	}
	else
	{
		ind->obj = NULL;
	}

    return;
}
