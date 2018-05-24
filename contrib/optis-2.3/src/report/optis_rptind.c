# include "../optis_global.h"

void report_ind(individual *ind, int nobj, FILE *file)
{
    int i;

    for(i=0;i<nobj;i++)
    {
    	fprintf(file, "%20e", ind->obj[i]);
    }

    fprintf(file, "%8d", ind->rank);
    fprintf(file, "%18e\n", ind->crowd_dist);

    return;
}
