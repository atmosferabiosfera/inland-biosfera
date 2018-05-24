# include <stdlib.h>

# include "../optis_global.h"

void dealloc_sitarray(sit_array *sitios)
{
    int i;

    for(i=0;i<sitios->size;i++)
    {
    	dealloc_sitio(&sitios->sit[i]);
    }

    free(sitios);

    return;
}
