# include <stdlib.h>

# include "optis_global.h"

/* O array obj do struct calibration_lvl possui 3 dimensoes, onde a primeira indica quais sitios precisam ser executados, *
 * a segunda indica em quais escalas temporais os dados serao analisados, e a terceira indica quais os metodos utilizados *
 * para calcular os valores das funcoes objetivos.									  */
void evaluate_ind(int ind_id, individual *ind, calibration_lvl *calib, base_param *bparam, sit_array *sitios, modelo_st *modelo)
{
	int i,k,l;
	double **output;
	double *of_value;
	var_array *var_list;

	for(i=0; i<calib->nobj; i++){
		ind->obj[i] = 0.0;
	}

	for(i=0;i<sitios->size;i++)
	{
		if(calib->obj[i] != NULL)
		{
			run_inland(ind->xreal, calib, bparam, ind_id, i, modelo);

			output = rd_inlandout(ind_id, i, sitios->sit[i].nlines, sitios->sit[i].ncol_out, sitios->sit[i].exline_out, modelo);

			for(k=0;k<N_OF_METHODS;k++)
			{
				if(calib->obj[i][k].size != 0)
				{
					var_list = &calib->obj[i][k];

					if(k == MAE)
					{
						of_value = calc_mae(sitios->sit[i].obs_data, output, sitios->sit[i].nlines, var_list);
					} else if(k == RMSE)
					{
						of_value = calc_rmse(sitios->sit[i].obs_data, output, sitios->sit[i].nlines, var_list);
					}
					for(l=0;l<var_list->size;l++)
					{
						ind->obj[var_list->var[l].pos_obj] = ind->obj[var_list->var[l].pos_obj] + of_value[l];
					}

					free(of_value);
				}
			}

			deallocate2d(output, sitios->sit[i].nlines);
		}
	}

	for(i=0; i<calib->nobj; i++){
		ind->obj[i] = ind->obj[i] / calib->obj_comp[i];
	}

	return;
}
