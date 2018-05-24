# include "../optis_global.h"

FILE* open_calibparamsf()
{
    FILE *file;
    char name[500];

    sprintf(name, "%sALL.calibratedparams.out", OUTPUTDIR);
    file = fopen(name, "w");

    fprintf(file, "# This file contains the values of the calbrated params for each calibration level.\n");

    return file;
}
