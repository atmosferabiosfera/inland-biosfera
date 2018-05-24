# include <stdlib.h>

# include "../optis_global.h"

void dealloc_baseparam(base_param *bparam)
{
	int i;

	for(i=0;i<bparam->nfiles;i++)
	{
		free(bparam->filename[i]);
	}
	free(bparam->filename);

	free(bparam);

    return;
}
