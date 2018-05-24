# include "../optis_global.h"

void report_calibparams(calibration_lvl *calib, FILE *file, individual *ind, char **parname)
{
    int i;

    fprintf(file, "\n# calibration level = %d, calibration name =  %s\n", calib->lvl, calib->name);
    fprintf(file, "# number of objectives = %d, number of calibrated parameters = %d\n", calib->nobj, calib->cpars);

    fprintf(file, "%1s%30s%30s", "#", "name", "value\n");

    for(i=0;i<calib->cpars;i++)
    {
    	fprintf(file, "%5d%30s%30.14le\n", i, parname[i], ind->xreal[i]);
    }

    fprintf(file, "\n");
    fflush(file);

    return;
}
