# include <stdlib.h>

# include "../optis_global.h"

void deallocate2d(double **matrix, int nlines)
{
    int i;

    for(i=0;i<nlines;i++)
    {
    	free(matrix[i]);
    }

    free(matrix);

    return;
}
