# include "../optis_global.h"

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void copy_ind (individual *ind1, individual *ind2, int nreal, int nobj)
{
    int i;
    ind2->rank = ind1->rank;
    ind2->crowd_dist = ind1->crowd_dist;
    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            ind2->xreal[i] = ind1->xreal[i];
        }
    }
    for (i=0; i<nobj; i++)
    {
        ind2->obj[i] = ind1->obj[i];
    }
    return;
}
