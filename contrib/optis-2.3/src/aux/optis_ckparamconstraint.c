# include <string.h>
# include <stdlib.h>

# include "../optis_global.h"

void check_paramconstraint(char *param_name, int parnum, constraint_array* constraint)
{
	char *rootname, *s_pft;
	int pft, i;

	rootname = (char*)malloc(100*sizeof(char));
	memset(rootname, '\0', sizeof(rootname));

	s_pft = (char*)malloc(10*sizeof(char));
	memset(s_pft, '\0', sizeof(s_pft));

	strncpy(rootname, param_name, 5);

	if(strcmpe(rootname, "aroot") == 0)
	{
		printf("\n Parameter %s can not be calibrated! (aroot = 1 - awood - aleaf).\n", param_name);
		exit(EXIT_FAILURE);
	}
	else if( (strcmpe(rootname, "awood") == 0) || (strcmpe(rootname, "aleaf") == 0) )
	{
		strncpy(s_pft, param_name+5, 6);
		pft = atoi(s_pft);

		if(constraint->size == 0)
		{
			i = -1;
		}
		else
		{
			for (i = 0; i < constraint->size; i++) {
				if(constraint->cons[i][0] == pft)
				{
					break;
				}
			}
		}

		if( (i == -1) || (i == constraint->size) )
		{
			add_constraint(constraint);
			if(i == -1)
			{
				i++;
			}
			constraint->cons[i][0] = pft;
		}

		if(strcmpe(rootname, "awood") == 0)
		{
			constraint->cons[i][1] = parnum;
		}
		else
		{
			constraint->cons[i][2] = parnum;
		}
	}
	free(rootname);
	free(s_pft);

	return;

}
