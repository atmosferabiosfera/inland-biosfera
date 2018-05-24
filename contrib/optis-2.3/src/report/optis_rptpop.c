# include "../optis_global.h"

void report_pop(FILE *file, population *pop, int gen, int popsize, int nobj)
{
	int i;

	fprintf(file, "# gen = %d\n", gen);

	for(i=0;i<popsize;i++)
	{
		report_ind(&(pop->ind[i]), nobj, file);
	}

	fflush(file);

    return;
}
