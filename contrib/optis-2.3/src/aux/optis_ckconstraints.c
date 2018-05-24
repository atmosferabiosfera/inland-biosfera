# include <stdlib.h>

# include "../optis_global.h"

void check_constraints(constraint_array *constraint, base_param *bparam)
{
	int i;
	char name[100];

	if(constraint->size > 0)
	{
		for (i = 0; i < constraint->size; ++i) {
			if( (constraint->cons[i][1] == -1) || (constraint->cons[i][2] == -1) )
			{
				printf("\n The parameters awood%d and aleaf%d must be calibrated together!\n", constraint->cons[i][0], constraint->cons[i][0]);
				exit(EXIT_FAILURE);
			}
			else
			{
				sprintf(name, "aroot%d", constraint->cons[i][0]);
//				constraint->cons[i][3] = find_param(name, bparam);
			}
		}
	}
	return;
}
