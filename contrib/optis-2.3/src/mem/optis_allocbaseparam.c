# include <stdlib.h>

# include "../optis_global.h"

base_param* alloc_baseparam(int nfiles)
{
	base_param *bparam;

	bparam = (base_param*)malloc(sizeof(base_param));

	bparam->nfiles = nfiles;

	bparam->filename = (char**)malloc(nfiles*sizeof(char*));

	return bparam;
}
