# include <string.h>
# include <stdlib.h>

# include "../optis_global.h"

void report_inlandparams(base_param *bparams)
{
    char *cmd;
    int j;

    cmd = (char*)malloc(500*sizeof(char));

    for(j=0;j<bparams->nfiles;j++)
    {
    	sprintf(cmd, "cp %s%s/%s %s;", TEMPDIR, AUXBASEDIR, bparams->filename[j], OUTPUTDIR);
    	system(cmd);
    }

    free(cmd);

    return;
}
