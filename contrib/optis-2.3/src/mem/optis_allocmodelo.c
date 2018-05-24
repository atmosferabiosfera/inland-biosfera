# include <stdlib.h>
# include <string.h>

# include "../optis_global.h"

modelo_st* alloc_modelo(){

	modelo_st* modelo;

	modelo = (modelo_st*)malloc(sizeof(modelo_st));

	return modelo;
}
