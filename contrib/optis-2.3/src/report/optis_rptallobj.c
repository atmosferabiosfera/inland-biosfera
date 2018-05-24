# include <stdlib.h>

# include "../optis_global.h"

void report_allobj(calibration_lvl *allobj, base_param *bparam, sit_array *sitios, modelo_st *modelo)
{
	individual *ind;
	FILE *file;
	char name[500];
	int i;

	ind = (individual*)malloc(sizeof(individual));

	ind->rank = 1;
	ind->crowd_dist = 0;

	alloc_ind(ind, 0, allobj->nobj);

	evaluate_ind(0, ind, allobj, bparam, sitios, modelo);

	sprintf(name, "%sALL.obj.out", OUTPUTDIR);
	file = fopen(name, "w");

	fprintf(file, "# This file contains all objective functions calculated with the final set of parameters.\n");
	fprintf(file, "\n#calibration name =  %s\n", allobj->name);
	fprintf(file, "# number of objectives = %d\n", allobj->nobj);

	for (i = 0; i < allobj->nobj; ++i)
	{
		fprintf(file, "%20s", allobj->of_name[i]);
	}
	fprintf(file, "%8s%18s\n", "rank", "c_distance");

	report_ind(ind, allobj->nobj, file);

	fflush(file);
	fclose(file);

	dealloc_ind(ind);
	free(ind);

	return;
}
