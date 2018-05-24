# include "../optis_global.h"

void set_paramsrange(int cal, int popsize, population *pop, int npars, paramsrange_opt *prange)
{
    int i, j;

    init_paramsrange(prange, cal, npars);

    for(i=0;i<popsize;i++)
    {
    	if(pop->ind[i].rank == 1)
    	{
    		for(j=0;j<npars;j++)
    		{
    			if(pop->ind[i].xreal[j] > prange->max_value[cal][j])
    			{
    				prange->max_value[cal][j] = pop->ind[i].xreal[j];
    			}

    			if(pop->ind[i].xreal[j] < prange->min_value[cal][j])
    			{
    				prange->min_value[cal][j] = pop->ind[i].xreal[j];
    			}
    		}
    	}
    }

    return;
}
