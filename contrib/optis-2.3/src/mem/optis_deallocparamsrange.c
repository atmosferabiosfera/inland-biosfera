# include <stdlib.h>

# include "../optis_global.h"

void dealloc_paramsrange(paramsrange_opt *prange)
{
	int i;

	for(i=0;i<prange->ncal;i++)
	{
		free(prange->min_value[i]);
		free(prange->max_value[i]);
	}
	free(prange->min_value);
	free(prange->max_value);

	free(prange);

    return;
}
