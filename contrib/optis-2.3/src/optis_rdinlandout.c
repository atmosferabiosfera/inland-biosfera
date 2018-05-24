# include <string.h>
# include <stdlib.h>

# include "optis_global.h"

double** rd_inlandout(int ind_id, int sitio_id, int nlines, int ncol_out, int exline_out, modelo_st *modelo)
{
	char file[500];
	double **output;

	if(strcmp(modelo->out_name, "NULL") != 0)
	{
		sprintf(file, "%s%d/%s%d/%s", TEMPDIR, ind_id, B_SITIODIR, sitio_id, modelo->out_name);
		output = rd_filewdouble(file, nlines, ncol_out, exline_out);
	}
	else
	{
		output = NULL;
		printf("ERRO = Arquivo de output %s não encontrado!\n", modelo->out_name);
		exit(EXIT_FAILURE);
	}

    return output;
}
