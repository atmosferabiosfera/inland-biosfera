/* Funcao para escrever os parametros do array aux_params em formato    *
 * dos arquivos de parametros padrao (bparams).				*
 * PATH e o local onde os arquivos serao escritos.			*/
# include <stdlib.h>

# include "optis_global.h"

void write_inlandparams(char *path, double *xreal, base_param *bparams, calibration_lvl *calib)
{
    int i, j;
    char cmd[600], *aux;
    int pf, pl, pc, ps;

    for(i=0;i<calib->cpars;i++)
    {
    	pf = calib->parf[i];
    	pl = calib->parl[i];
    	pc = calib->parc[i];
    	ps = calib->pars[i];

    	aux = (char*)malloc((ps+1)*sizeof(char));
    	for(j=0;j<ps;j++)
    	{
    		aux[j] = '.';
    	}
    	aux[ps] = '\0';

    	sprintf(cmd, "mv %s%s %spar_temp; sed '%ds/^\\(.\\{%d\\}\\)%s\\(.*\\)/\\1%+.*E\\2/' < %spar_temp > %s%s;rm %spar_temp;",
    			path, bparams->filename[pf], path, pl, pc-1, aux, ps-7, xreal[i], path, path, bparams->filename[pf], path);

    	system(cmd);

    	free(aux);
    }

    return;
}
