# include <stdlib.h>

# include "../optis_global.h"

calibration_lvl* alloc_allobj(sit_array *sitios)
{
	calibration_lvl *allobj;
	int j;

	allobj = (calibration_lvl*)malloc(sizeof(calibration_lvl));

	allobj->cpars = 0;
	allobj->parl = NULL;
	allobj->parc = NULL;
	allobj->pars = NULL;
	allobj->parf = NULL;
	allobj->parname = NULL;
	allobj->min_parvalue = NULL;
	allobj->max_parvalue = NULL;
	allobj->nobj = 0;
	allobj->of_name = NULL;
	allobj->lvl = 0;
	allobj->obj_comp = NULL;

	allobj->constraint = alloc_constraint();

	allobj->obj = (var_array**)malloc(sitios->size*sizeof(var_array*));
    for (j = 0; j < sitios->size; j++) {
		allobj->obj[j] = NULL;
	}

	allobj->name = (char*)malloc(100*sizeof(char));
	sprintf(allobj->name, "All objectives");

	return allobj;
}
