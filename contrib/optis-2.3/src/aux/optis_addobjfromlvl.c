# include <string.h>
# include <stdlib.h>

# include "../optis_global.h"

void add_objfromlvl(calibration_lvl *allobj, calibration_lvl *calib, sit_array *sitios)
{
	int i, j, sit, method;
	var_array *var_list;
	variavel* var;

	for(sit=0;sit<sitios->size;sit++)
	{
		if(calib->obj[sit] != NULL)
		{
			if(allobj->obj[sit] == NULL)
			{
				allobj->obj[sit] = (var_array*) malloc(N_OF_METHODS * sizeof(var_array));
				for (j = 0; j < N_OF_METHODS; j++)
				{
					allobj->obj[sit][j].size = 0;
				}
			}
			for(method=0;method<N_OF_METHODS;method++)
			{
				if(calib->obj[sit][method].size != 0)
				{
					var_list = &calib->obj[sit][method];
					for(i=0;i<var_list->size;i++)
					{
						var = (variavel*) malloc(sizeof(variavel));
						var->pos_obs = var_list->var[i].pos_obs;
						var->pos_out = var_list->var[i].pos_out;
						var->pos_obj = allobj->nobj + var_list->var[i].pos_obj;

						add_variavel(&allobj->obj[sit][method], var);

					}
				}
			}
		}
	}

	if(allobj->nobj == 0)
	{
		allobj->of_name = (char**)malloc(calib->nobj * sizeof(char*));
		allobj->obj_comp = (int*)malloc(calib->nobj * sizeof(int));
	}
	else
	{
		allobj->of_name = (char**)realloc(allobj->of_name, (allobj->nobj + calib->nobj)*sizeof(char*));
		allobj->obj_comp = (int*)realloc(allobj->obj_comp, (allobj->nobj + calib->nobj)*sizeof(int));
	}

	for(i=0; i < calib->nobj; i++){
		allobj->of_name[allobj->nobj + i] = (char*)malloc(100*sizeof(char));
		strcpy(allobj->of_name[allobj->nobj + i], calib->of_name[i]);
		allobj->obj_comp[allobj->nobj + i] = calib->obj_comp[i];
	}

	allobj->nobj = allobj->nobj + calib->nobj;

	return;
}
