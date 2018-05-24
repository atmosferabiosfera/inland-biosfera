# include <stdlib.h>
# include <string.h>

# include "../optis_global.h"

void rd_calibrationvars(int lvl, calibration_lvl* calib, FILE* file, sit_array* sitios) {

	int i, sit, j, method;
	char method_name[100], *var_name;
	variavel* var;
	int aux_nobj, *aux_obj_comp, pos_obj;
	char **aux_name, **aux_of_name;

	calib->nobj = 0;
	aux_nobj = rd_int(file);
	if(aux_nobj < 1){
		printf("\nERRO: O total de funcoes objetivas precisa ser maior que zero.\n");
		exit(EXIT_FAILURE);
	}

	aux_name = (char**) malloc(aux_nobj * sizeof(char*));
	aux_obj_comp = (int*) malloc(aux_nobj * sizeof(int));
	aux_of_name = (char**) malloc(aux_nobj * sizeof(char*));

	for(i = 0; i < aux_nobj; i++){
		aux_obj_comp[i] = 0;
		aux_of_name[i] = NULL;
	}

	for (i = 0; i < aux_nobj; i++) {

		var_name = rd_string(file, 100);

		if(calib->mode == CALIB_MOD_MED){
			//Buscando por variavel de mesmo nome
			for(j = 0; j < calib->nobj; j++){
				if(strcmp(aux_name[j], var_name) == 0){
					pos_obj = j;
					break;
				}
			}
		}

		//Variavel nova
		if( (calib->mode == CALIB_MOD_IND) || (j == calib->nobj) ){
			aux_name[calib->nobj] = (char*) malloc(100 * sizeof(char));
			strcpy(aux_name[calib->nobj], var_name);
			pos_obj = calib->nobj;
			calib->nobj = calib->nobj + 1;
		}

		aux_obj_comp[pos_obj] = aux_obj_comp[pos_obj] + 1;

		/* Lendo o sitio ao qual a variavel pertence */
		sit = rd_int(file) - 1;
		if (sit >= sitios->size || sit < 0) {
			printf("\nERRO: O sitio %d da variavel %s no nivel %d nao existe.\n", sit+1, var_name, lvl);
			exit(EXIT_FAILURE);
		} else {
			if (calib->obj[sit] == NULL) {
				calib->obj[sit] = (var_array*) malloc(N_OF_METHODS * sizeof(var_array));
				for (j = 0; j < N_OF_METHODS; j++) {
					calib->obj[sit][j].size = 0;
				}
			}
		}


		/* Lendo o metodo para calculo da funcao objetiva */
		method = rd_int(file);
		if (method >= N_OF_METHODS || method < 0) {
			printf("\nERRO: O metodo %d da variavel %s no nivel %d nao existe.\n", method, var_name, lvl);
			exit(EXIT_FAILURE);
		}

		var = (variavel*) malloc(sizeof(variavel));

		var->pos_obs = rd_int(file) - 1;
		if (var->pos_obs > sitios->sit[sit].ncol_obs || var->pos_obs < 1) {
			printf("\nERRO: A posicao da variavel %s (nivel %d) no arquivo de dados observado precisa estar entre 1 e %d.\n", var_name, lvl, sitios->sit[sit].ncol_obs);
			exit(EXIT_FAILURE);
		}

		var->pos_out = rd_int(file) - 1;
		if (var->pos_out > sitios->sit[sit].ncol_out || var->pos_out < 1) {
			printf("\nERRO: A posicao da variavel %s (nivel %d) no arquivo de output precisa estar entre 1 e %d.\n", var_name, lvl, sitios->sit[sit].ncol_out);
			exit(EXIT_FAILURE);
		}

		var->pos_obj = pos_obj;

		/* Add a variavel ao conjunto de variaves */
		add_variavel(&calib->obj[sit][method], var);

		/* Testar com todos os metodos existentes.*/
		if (method == MAE) {
			strcpy(method_name, MAE_STRING);
		}else if (method == RMSE) {
			strcpy(method_name, RMSE_STRING);
		}


		if(aux_of_name[pos_obj] == NULL){

			aux_of_name[pos_obj] = (char*) malloc(100 * sizeof(char));

			if(calib->mode == CALIB_MOD_IND){
			    sprintf(aux_of_name[pos_obj], "%s(%s):%s", method_name, var_name, sitios->sit[sit].name);
			}else{
				sprintf(aux_of_name[pos_obj], "%s(%s):m", method_name, var_name);
			}
		}

		free(var_name);

	}

	calib->of_name = (char**) malloc(calib->nobj * sizeof(char*));
	calib->obj_comp = (int*) malloc(calib->nobj * sizeof(int));

	for(i = 0; i < calib->nobj; i++){

		calib->of_name[i] = (char*) malloc(100 * sizeof(char));
		strcpy(calib->of_name[i], aux_of_name[i]);

		calib->obj_comp[i] = aux_obj_comp[i];

		free(aux_name[i]);
		free(aux_of_name[i]);

	}

	free(aux_obj_comp);
	free(aux_name);
	free(aux_of_name);
}
