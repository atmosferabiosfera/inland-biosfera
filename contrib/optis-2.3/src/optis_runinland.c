# include <stdlib.h>
# include <string.h>

# include "optis_global.h"

void run_inland(double *xreal, calibration_lvl *calib, base_param *bparams, int ind_id, int sitio_id, modelo_st *modelo)
{
	char *path, cmd[300];
	constraint_array *constraint;

	constraint = calib->constraint;

//	if(constraint->size > 0)
//	{
//		for (i = 0; i < constraint->size; ++i) {
//			aux_params[constraint->cons[i][3]] = 1.0 - (aux_params[constraint->cons[i][1]] + aux_params[constraint->cons[i][2]]);
//		}
//	}

	path = (char*)malloc(500*sizeof(char));

	sprintf(path, "%s%d/%s%d/", TEMPDIR, ind_id, B_SITIODIR, sitio_id);

	write_inlandparams(path, xreal, bparams, calib);

	sprintf(cmd, "cd %s%d/%s%d/; ./%s >> /dev/null", TEMPDIR, ind_id, B_SITIODIR, sitio_id, modelo->bin);

	system(cmd);

	free(path);

	return;

}
