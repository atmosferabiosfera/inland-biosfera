# include <stdlib.h>

# include "../optis_global.h"

nsga_param* init_nsgaparam(FILE *file)
{
    nsga_param *nparam;

    nparam = alloc_nsgaparam();

    nparam->ngen = rd_int(file);
    nparam->popsize = rd_int(file);

    if((nparam->ngen < 1) || (nparam->popsize < 1)){
    	printf("\nERRO: O numero de geracoes e tamanho da populacao deve ser maior que zero.\n");
    	exit(EXIT_FAILURE);
    }

    nparam->pcross_real = rd_double(file);
    nparam->pmut_real = rd_double(file);

    if((nparam->pcross_real < 0) || (nparam->pcross_real > 1) || (nparam->pmut_real < 0) || (nparam->pmut_real > 1)){
    	printf("\nERRO: pcross_real e pmut_real devem estar entre 0 e 1.\n");
    	exit(EXIT_FAILURE);
    }

    nparam->eta_c = rd_double(file);
    nparam->eta_m = rd_double(file);
    nparam->seed = rd_double(file);

    if((nparam->eta_c < 0) || (nparam->eta_m < 0) || (nparam->seed < 0)){
    	printf("\nERRO: eat_c, eta_m e seed devem ser positivos.\n");
    	exit(EXIT_FAILURE);
    }

    return nparam;
}
