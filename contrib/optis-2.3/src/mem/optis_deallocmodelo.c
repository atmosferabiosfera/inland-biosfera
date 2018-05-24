# include <stdlib.h>

# include "../optis_global.h"

void dealloc_modelo(modelo_st *modelo){

	free(modelo->bin);
	free(modelo->st_script);
	free(modelo->out_name);

	free(modelo);

	return;

}
