/* Mutation routines */

# include <math.h>

# include "../optis_global.h"

/* Function to perform mutation in a population */
void mutation_pop (population *pop, int popsize, int nreal, double pmut_real, double *min_realvar, double *max_realvar, double eta_m)
{
    int i, j;

    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    individual *ind;

    for (i=0; i<popsize; i++)
    {
	ind = &(pop->ind[i]);

	for (j=0; j<nreal; j++)
	{
	    if (randomperc() <= pmut_real)
	    {
		y = ind->xreal[j];
		yl = min_realvar[j];
		yu = max_realvar[j];
		delta1 = (y-yl)/(yu-yl);
		delta2 = (yu-y)/(yu-yl);
		rnd = randomperc();
		mut_pow = 1.0/(eta_m+1.0);
		if (rnd <= 0.5)
		{
		    xy = 1.0-delta1;
		    val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
		    deltaq =  pow(val,mut_pow) - 1.0;
		}
		else
		{
		    xy = 1.0-delta2;
		    val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
		    deltaq = 1.0 - (pow(val,mut_pow));
		}
		y = y + deltaq*(yu-yl);
		if (y<yl)
		    y = yl;
		if (y>yu)
		    y = yu;
		ind->xreal[j] = y;
		/*nrealmut+=1;*/
	    }
	}
    }
    return;
}

