# include <stdlib.h>

# include "../optis_global.h"

constraint_array* alloc_constraint()
{
	constraint_array *constraint;

	constraint = (constraint_array*)malloc(sizeof(constraint_array));

	constraint->size = 0;

	constraint->cons = NULL;

	return constraint;
}
