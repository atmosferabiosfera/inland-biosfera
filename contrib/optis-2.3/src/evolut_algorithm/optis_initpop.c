/* Data initializtion routines */

# include "../optis_global.h"

/* Function to initialize a population randomly */
void initialize_pop (population *pop, int popsize, int nreal, double *min_realvar, double *max_realvar)
{
    int i, j;

    for (i=0; i<popsize; i++)
    {
        for (j=0; j<nreal; j++)
        {
            pop->ind[i].xreal[j] = rndreal (min_realvar[j], max_realvar[j]);
        }

    }

    return;
}

