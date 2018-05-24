/* Routine for mergeing two populations */

# include "../optis_global.h"

/* Routine to merge two populations into one */
void merge(population *pop1, population *pop2, population *pop3, int popsize, int nreal, int nobj)
{
    int i, k;
    for (i=0; i<popsize; i++)
    {
        copy_ind (&(pop1->ind[i]), &(pop3->ind[i]), nreal, nobj);
    }
    for (i=0, k=popsize; i<popsize; i++, k++)
    {
        copy_ind (&(pop2->ind[i]), &(pop3->ind[k]), nreal, nobj);
    }
    return;
}
