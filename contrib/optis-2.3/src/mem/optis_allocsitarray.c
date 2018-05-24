# include <stdlib.h>

# include "../optis_global.h"

sit_array* alloc_sitarray(int size)
{
    int i;
    sit_array *sitios;

    sitios = (sit_array*)malloc(sizeof(sit_array));

    sitios->size = size;

    sitios->sit = (sitio*)malloc(size*sizeof(sitio));

    for(i=0;i<size;i++)
    {
    	alloc_sitio(&(sitios->sit[i]));
    }

    return sitios;
}
