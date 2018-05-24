# include <math.h>
# include <stdlib.h>

# include "../optis_global.h"

double* calc_rmse(double **obs, double **output, int nlines, var_array *var_list)
{
	double *rmse;
    int i, j, n_vars;
    int *count;

    n_vars = var_list->size;

    rmse = (double*)malloc(n_vars*sizeof(double));
    count = (int*)malloc(n_vars*sizeof(int));

    for(i=0;i<n_vars;i++)
    {
    	rmse[i] = 0.0;
    	count[i] = 0;
    }

    for(i=0;i<nlines;i++)
    {
    	for(j=0;j<n_vars;j++)
    	{
    		if( obs[i][var_list->var[j].pos_obs] != ND_VALUE )
    		{
    			rmse[j] = rmse[j] + pow(obs[i][var_list->var[j].pos_obs] - output[i][var_list->var[j].pos_out], 2.0);
    			count[j]++;
    		}
    	}
    }

    for(i=0;i<n_vars;i++)
    {
    	if(count[i] == 0)
    	{
    		rmse[i] = ND_VALUE;
    	}
    	else
    	{
    		rmse[i] = sqrt(rmse[i] / count[i]);
    	}
    }

    free(count);

    return rmse;



}
