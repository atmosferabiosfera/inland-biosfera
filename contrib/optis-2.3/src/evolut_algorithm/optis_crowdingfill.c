# include <stdlib.h>

# include "../optis_global.h"

/* Routine to fill a population with individuals in the decreasing order of crowding distance */
void crowding_fill (population *mixed_pop, population *new_pop, int count, int front_size, list *elite, int nobj, int popsize, int nreal)
{
    int *dist;
    list *temp;
    int i, j;
    assign_crowding_distance_list (mixed_pop, elite->child, front_size, nobj);
    dist = (int *)malloc(front_size*sizeof(int));
    temp = elite->child;
    for (j=0; j<front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }
    quicksort_dist (mixed_pop, dist, 0, front_size-1);
    for (i=count, j=front_size-1; i<popsize; i++, j--)
    {
        copy_ind(&mixed_pop->ind[dist[j]], &new_pop->ind[i], nreal, nobj);
    }
    free (dist);
    return;
}
