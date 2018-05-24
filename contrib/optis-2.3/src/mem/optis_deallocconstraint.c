# include <stdlib.h>

# include "../optis_global.h"

void dealloc_constraint(constraint_array *constraint)
{
	int i;

	for (i = 0; i < constraint->size; ++i) {
		free(constraint->cons[i]);
	}

	free(constraint->cons);

	free(constraint);

	return;
}
