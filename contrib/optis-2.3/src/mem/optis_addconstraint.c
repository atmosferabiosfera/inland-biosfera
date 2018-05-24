# include <stdlib.h>

# include "../optis_global.h"

void add_constraint(constraint_array *constraint)
{
	int i;

	if(constraint->size == 0)
	{
		constraint->cons = (int**)malloc(sizeof(int*));
		constraint->cons[0] = (int*)malloc(4*sizeof(int));

		for (i = 0; i < 4; ++i) {
			constraint->cons[0][i] = -1;
		}

		constraint->size = 1;
	}
	else
	{
		constraint->size++;
		constraint->cons = (int**)realloc(constraint->cons, constraint->size*sizeof(int*));

		constraint->cons[constraint->size - 1] = (int*)malloc(4*sizeof(int));

		for (i = 0; i < 4; ++i) {
			constraint->cons[constraint->size - 1][i] = -1;
		}

	}
	return;
}
