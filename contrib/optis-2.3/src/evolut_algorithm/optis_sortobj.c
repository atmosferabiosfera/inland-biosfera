# include "../optis_global.h"

/* Randomized quick sort routine to sort a population based on a particular objective chosen */
void quicksort_front_obj(population *pop, int objcount, int obj_array[], int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;

    if (left<right)
    {
        /*index = rnd (left, right);*/
    	index = (left + right) / 2;
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        pivot = pop->ind[obj_array[right]].obj[objcount];
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (pop->ind[obj_array[j]].obj[objcount] <= pivot)
            {
                i+=1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }
        index=i+1;
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        quicksort_front_obj(pop, objcount, obj_array, left, index-1);
        quicksort_front_obj(pop, objcount, obj_array, index+1, right);
    }

    return;
}


