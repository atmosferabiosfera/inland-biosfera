# include <stdlib.h>

# include "../optis_global.h"

void dealloc_pop(population *pop, int popsize)
{
	int i;

	for(i=0;i<popsize;i++)
	{
		dealloc_ind(&(pop->ind[i]));
	}

	free(pop->ind);

	free(pop);

    return;
}
