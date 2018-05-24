# include "../optis_global.h"

void init_paramsrange(paramsrange_opt *prange, int cal, int npars)
{
	int j;

	for(j=0;j<npars;j++)
	{
		prange->min_value[cal][j] = 999999;
		prange->max_value[cal][j] = -999999;
	}

	return;
}
