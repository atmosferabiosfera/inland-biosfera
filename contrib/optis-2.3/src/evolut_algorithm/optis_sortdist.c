# include "../optis_global.h"

/* Randomized quick sort routine to sort a population based on crowding distance */
void quicksort_dist(population *pop, int *dist, int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;

    if (left<right)
    {
        /*index = rnd (left, right);*/
        index = (left + right) / 2;
        temp = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = pop->ind[dist[right]].crowd_dist;
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (pop->ind[dist[j]].crowd_dist <= pivot)
            {
                i+=1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            }
        }
        index=i+1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        quicksort_dist(pop, dist, left, index-1);
        quicksort_dist(pop, dist, index+1, right);
    }

    return;
}


