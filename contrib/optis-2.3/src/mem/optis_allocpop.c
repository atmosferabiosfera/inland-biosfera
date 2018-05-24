# include <stdlib.h>

# include "../optis_global.h"

population* alloc_pop(int popsize, int cpars, int nobj)
{
    int i;
    population *pop;

    pop = (population*)malloc(sizeof(population));

    pop->ind = (individual*)malloc(popsize*sizeof(individual));

    for(i=0;i<popsize;i++)
    {
    	alloc_ind(&(pop->ind[i]), cpars, nobj);
    }

    return pop;
}
