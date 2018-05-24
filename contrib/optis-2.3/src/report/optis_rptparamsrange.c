# include "../optis_global.h"

void report_paramsrange(paramsrange_opt *prange, calibration_lvl **allcaliblvl)
{
    int i, j;
    FILE *file;
    char name[500];
    calibration_lvl *calib;

    sprintf(name, "%sALL.paramsrange.out", OUTPUTDIR);
    file = fopen(name, "w");

    fprintf(file, "# This file contains the range of the calibrated parameters from non-dominated individuals.\n");

    for(i=0;i<prange->ncal;i++)
    {
    	calib = allcaliblvl[i];

    	fprintf(file, "\n# calibration level = %d, calibration name =  %s\n", calib->lvl, calib->name);
    	fprintf(file, "# number of calibrated parameters = %d\n", calib->cpars);
    	fprintf(file, "%30s%30s%30s", "param", "min", "max\n");

    	for(j=0;j<calib->cpars;j++)
    	{
    		fprintf(file, "%30s%30.14le%30.14le\n", calib->parname[j], prange->min_value[i][j], prange->max_value[i][j]);
    	}
    }

    fprintf(file, "\n");
    fflush(file);
    fclose(file);

    return;
}
