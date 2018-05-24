# include <stdlib.h>

# include "../optis_global.h"

void dealloc_allcalibrationlvl(int ncal, calibration_lvl **allcaliblvl, sit_array *sitios)
{
	int i;

	for(i=0;i<ncal;i++)
	{
		dealloc_calibrationlvl(allcaliblvl[i], sitios);
	}

	free(allcaliblvl);

    return;
}
