# include <stdlib.h>

# include "../optis_global.h"

/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */
void assign_crowding_distance_list (population *pop, list *lst, int front_size, int nobj)
{
    int **obj_array;
    int *dist;
    int i, j;
    list *temp;
    temp = lst;
    if (front_size==1)
    {
        pop->ind[lst->index].crowd_dist = INF;
        return;
    }
    if (front_size==2)
    {
        pop->ind[lst->index].crowd_dist = INF;
        pop->ind[lst->child->index].crowd_dist = INF;
        return;
    }
    obj_array = (int **)malloc(nobj*sizeof(int*));
    dist = (int *)malloc(front_size*sizeof(int));
    for (i=0; i<nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size*sizeof(int));
    }
    for (j=0; j<front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }
    assign_crowding_distance (pop, dist, obj_array, front_size, nobj);
    free (dist);
    for (i=0; i<nobj; i++)
    {
        free (obj_array[i]);
    }
    free (obj_array);
    return;
}
