/*
 * optis_rptparetofront.c
 *
 *  Created on: 23/02/2017
 *      Author: vitor
 */

# include <stdlib.h>

# include "../optis_global.h"

void report_paretofront(base_param *bparam, int popsize, population *pop, calibration_lvl *calib)
{
	char cmd[500], dir[200], name[200], inddir[200], path[500];
	FILE *file;
	int i, j;

	/* Cria, se ainda nao existe, o diretorio dos dados da fronteira de pareto */
	sprintf(cmd, "mkdir -p %s;", OUTPARETODIR);
	system(cmd);

	/* Cria o diretorio para o nivel hierarquico atual */
    sprintf(dir, "%snivel%d/", OUTPARETODIR, calib->lvl);
    sprintf(cmd, "mkdir %s;", dir);
	system(cmd);


	sprintf(name, "%sfront.list", dir);
	file = fopen(name, "w");

	fprintf(file, "# This file contains the list of individuals in the Pareto Front.\n\n");

	fprintf(file, "%5s", "ind");
	for(i=0;i<calib->nobj;i++)
	{
		fprintf(file, "%20s", calib->of_name[i]);
	}
	fprintf(file, "%10s\n", "rank");


	for(i=0;i<popsize;i++)
	{
		if(pop->ind[i].rank == 1)
	   	{
			fprintf(file, "%5d", i);

			for(j=0;j<calib->nobj;j++)
			{
			   	fprintf(file, "%20e", pop->ind[i].obj[j]);
			}

			fprintf(file, "%10d\n", pop->ind[i].rank);

			/* Cria o diretorio para o individuo da fronteira */
			sprintf(inddir, "%sind_%d/", dir, i);
			sprintf(cmd, "mkdir %s;", inddir);
		    system(cmd);

//	    	if(calib->constraint->size > 0)
//		    {
//		    	for (j = 0; j < calib->constraint->size; ++j) {
//		    		aux_params[calib->constraint->cons[j][3]] = 1.0 - (aux_params[calib->constraint->cons[j][1]] + aux_params[calib->constraint->cons[j][2]]);
//		    	}
//		    }

		    /* Imprime os arquivos de parametros do individuo */
		    sprintf(path, "%s%s/", TEMPDIR, AUXBASEDIR);
		    write_inlandparams(path, pop->ind[i].xreal, bparam, calib);

		    for(j=0;j<bparam->nfiles;j++)
		    {
		    	sprintf(cmd, "cp %s%s/%s %s;", TEMPDIR, AUXBASEDIR, bparam->filename[j], inddir);
		    	system(cmd);
		    }

    	}
	}



	fflush(file);
	fclose(file);

	return;

}
