# include <stdlib.h>
# include <unistd.h>

# include "../optis_global.h"

base_param* init_baseparam(FILE *file)
{
    base_param *bparam;
    int nfiles, i;

    nfiles = rd_int(file);

    if(nfiles < 1){
    	printf("\nERRO: O total de arquivos de parametros do modelo precisa ser maior que zero.\n");
    	exit(EXIT_FAILURE);
    }

    bparam = alloc_baseparam(nfiles);

    for(i=0;i<nfiles;i++)
    {
    	bparam->filename[i] = rd_string(file, 100);
    }

    return bparam;
}
