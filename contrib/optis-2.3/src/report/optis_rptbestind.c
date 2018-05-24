# include "../optis_global.h"

void report_bestind(FILE *file, individual *ind, calibration_lvl *calib)
{
	int i;

	fprintf(file, "\n# calibration level = %d, calibration name = %s.\n", calib->lvl, calib->name);
	fprintf(file, "# number of objectives = %d, number of parameters = %d.\n", calib->nobj, calib->cpars);

	for(i=0;i<calib->nobj;i++)
	{
		fprintf(file, "%20s", calib->of_name[i]);
	}
	fprintf(file, "%8s%18s\n", "rank", "c_distance");

	report_ind(ind, calib->nobj, file);

	fflush(file);

	return;
}
