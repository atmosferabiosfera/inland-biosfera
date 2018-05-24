# include <stdlib.h>

# include "optis_global.h"

double** rd_filewdouble(char *name, int nline, int ncol, int exline)
{
	FILE *file;
	double **data;
	int i, j;

	file = fopen(name, "r");
	if(file == NULL)
	{
		return NULL;
	}

	if(exline > 0)
	{
		for(i=0;i<exline;i++)
		{
			for(j=0;j<ncol;j++)
			{
				fscanf(file, "%*s");
			}
		}
	}

	data = (double**)malloc(nline*sizeof(double*));
	for(i=0;i<nline;i++)
	{
		data[i] = (double*)malloc(ncol*sizeof(double));

		for(j=0;j<ncol;j++)
		{
			fscanf(file, "%lf", &data[i][j]);
		}
	}

	fclose(file);

	return data;

}
