# include <math.h>

# include "optis_global.h"

double calcdist_origin(individual *ind, double *maxvalue, int ndims)
{
	int i;
	double dist;

	dist = 0.0;
	for(i=0;i<ndims;i++)
	{
		if(maxvalue[i] > 0.01){
		    dist = dist + ( (ind->obj[i]/maxvalue[i]) * (ind->obj[i]/maxvalue[i]) );
		}
	}

	dist = sqrt(dist);

	return dist;
}
