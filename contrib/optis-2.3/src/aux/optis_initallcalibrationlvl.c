# include <stdlib.h>

# include "../optis_global.h"

calibration_lvl** init_allcalibrationlvl(FILE *f_conf, base_param *bparam, sit_array *sitios, int ncal, calibration_lvl *allobj)
{
	calibration_lvl **allcaliblvl;
	int i;

	allcaliblvl = (calibration_lvl**)malloc(ncal*sizeof(calibration_lvl*));

	for(i=0;i<ncal;i++)
	{
		allcaliblvl[i] = init_calibrationlvl(f_conf, i+1, bparam, sitios);

		add_objfromlvl(allobj, allcaliblvl[i], sitios);
	}

    return allcaliblvl;
}
