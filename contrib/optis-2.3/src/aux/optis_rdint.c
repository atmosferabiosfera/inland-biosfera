# include <stdlib.h>

# include "../optis_global.h"

int rd_int(FILE *file)
{
    int valor;
    char *aux;

    aux = rd_string(file, 100);

    valor = atoi(aux);

    free(aux);

    return valor;
}
