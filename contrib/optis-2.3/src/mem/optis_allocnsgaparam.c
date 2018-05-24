# include <stdlib.h>

# include "../optis_global.h"

nsga_param* alloc_nsgaparam()
{
    nsga_param *nparam;

    nparam = (nsga_param*)malloc(sizeof(nsga_param));

    return nparam;
}
