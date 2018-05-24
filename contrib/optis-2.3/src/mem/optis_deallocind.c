# include <stdlib.h>

# include "../optis_global.h"

void dealloc_ind(individual *ind)
{
    free(ind->xreal);

    free(ind->obj);

    return;
}
