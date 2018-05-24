# include <stdlib.h>

# include "../optis_global.h"

/* Routine for tournament selection, it creates a new_pop from old_pop by performing tournament selection and the crossover */
void selection (population *old_pop, population *new_pop, int popsize, int nobj, int nreal, double pcross_real, double *min_realvar, double *max_realvar, double eta_c)
{
    int *a1, *a2;
    int temp;
    int i;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *)malloc(popsize*sizeof(int));
    a2 = (int *)malloc(popsize*sizeof(int));
    for (i=0; i<popsize; i++)
    {
        a1[i] = a2[i] = i;
    }
    for (i=0; i<popsize; i++)
    {
        rand = rnd (i, popsize-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rnd (i, popsize-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    for (i=0; i<popsize; i+=4)
    {
        parent1 = tournament (&old_pop->ind[a1[i]], &old_pop->ind[a1[i+1]], nobj);
        parent2 = tournament (&old_pop->ind[a1[i+2]], &old_pop->ind[a1[i+3]], nobj);
        crossover (parent1, parent2, &new_pop->ind[i], &new_pop->ind[i+1], nreal, pcross_real, min_realvar, max_realvar, eta_c);
        parent1 = tournament (&old_pop->ind[a2[i]], &old_pop->ind[a2[i+1]], nobj);
        parent2 = tournament (&old_pop->ind[a2[i+2]], &old_pop->ind[a2[i+3]], nobj);
        crossover (parent1, parent2, &new_pop->ind[i+2], &new_pop->ind[i+3], nreal, pcross_real, min_realvar, max_realvar, eta_c);
    }
    free (a1);
    free (a2);
    return;
}
