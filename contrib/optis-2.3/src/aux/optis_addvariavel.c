# include <stdlib.h>

# include "../optis_global.h"

void add_variavel(var_array *varlist, variavel *var)
{
	if(varlist->size == 0)
	{
		varlist->var = var;
		varlist->size++;
	}
	else
	{
		varlist->size++;
		varlist->var = (variavel*)realloc(varlist->var, varlist->size*sizeof(variavel));

		varlist->var[varlist->size-1].pos_obs = var->pos_obs;
		varlist->var[varlist->size-1].pos_out = var->pos_out;
		varlist->var[varlist->size-1].pos_obj = var->pos_obj;

		free(var);
	}

    return;
}
