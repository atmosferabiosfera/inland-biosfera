# include <stdlib.h>

# include "../optis_global.h"

void dealloc_sitio(sitio *sit)
{
    int i, j;

    for(i=0;i<sit->nop;i++)
    {
    	for(j=0;j<2;j++)
    	{
    		free(sit->op[i][j]);
    	}
    	free(sit->op[i]);
    }
    free(sit->op);

    if(sit->obs_data != NULL){
   		deallocate2d(sit->obs_data, sit->nlines);
   	}

    free(sit->name);

    return;
}
