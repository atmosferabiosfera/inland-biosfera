# include <stdlib.h>

# include "../optis_global.h"

int** rd_filewint(char *name, int nline, int ncol, int exline)
{
	FILE *file;
	int **data;
	int i,j;

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

	data = (int**)malloc(nline*sizeof(int*));
	for(i=0;i<nline;i++)
	{
		data[i] = (int*)malloc(ncol*sizeof(int));

		for(j=0;j<ncol;j++)
		{
			fscanf(file, "%d", &data[i][j]);
		}
	}

	fclose(file);

	return data;
}
