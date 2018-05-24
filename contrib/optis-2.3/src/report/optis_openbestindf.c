# include "../optis_global.h"

FILE* open_bestindf()
{
    FILE *file;
    char name[500];

    sprintf(name, "%sALL.best_ind.out", OUTPUTDIR);
    file = fopen(name, "w");

    fprintf(file, "#This file contains the selected (optimal) individuals, from each calibration level.\n");

    return file;
}
