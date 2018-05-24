# include <stdlib.h>

# include "../optis_global.h"

void dealloc_nsgaparam(nsga_param* nparam)
{
    free(nparam);

    return;
}
