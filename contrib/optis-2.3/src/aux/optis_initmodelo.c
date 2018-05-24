# include <stdlib.h>
# include <unistd.h>

# include "../optis_global.h"

modelo_st* init_modelo(FILE *file){

	modelo_st* modelo;

	modelo = alloc_modelo();

	modelo->st_script = rd_string(file, 300);

	if(access(modelo->st_script, X_OK) != 0){
		printf("\nERRO: O arquivo %s nao existe ou nao pode ser executado.\n", modelo->st_script);
		exit(EXIT_FAILURE);
	}

	modelo->bin = rd_string(file, 100);

	/* Lendo o nome do arquivo de output do Modelo calibrado*/
	modelo->out_name = rd_string(file, 100);

	return modelo;

}
