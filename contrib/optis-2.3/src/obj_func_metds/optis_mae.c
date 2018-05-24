# include <math.h>
# include <stdlib.h>

# include "../optis_global.h"

double* calc_mae(double **obs, double **output, int nlines, var_array *var_list)
{
    double *mae;
    int i, j, n_vars;
    int *count;

    n_vars = var_list->size;

    mae = (double*)malloc(n_vars*sizeof(double));
    count = (int*)malloc(n_vars*sizeof(int));

    for(i=0;i<n_vars;i++)
    {
    	mae[i] = 0.0;
    	count[i] = 0;
    }

    for(i=0;i<nlines;i++)
    {
    	for(j=0;j<n_vars;j++)
    	{
    		if( obs[i][var_list->var[j].pos_obs] != ND_VALUE )
    		{
    			mae[j] = mae[j] + fabs(obs[i][var_list->var[j].pos_obs] - output[i][var_list->var[j].pos_out]);
    			count[j]++;
    		}
    	}
    }

    for(i=0;i<n_vars;i++)
    {
    	if(count[i] == 0)
    	{
    		mae[i] = ND_VALUE;
    	}
    	else
    	{
    		mae[i] = mae[i] / count[i];
    	}
    }

    free(count);

    return mae;
}
