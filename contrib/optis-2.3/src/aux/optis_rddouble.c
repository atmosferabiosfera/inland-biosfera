# include <stdlib.h>

# include "../optis_global.h"

double rd_double(FILE *file)
{
    double valor;
    char *aux;
    
    aux = rd_string(file, 100);

    valor = atof(aux);

    free(aux);

    return valor;

}
