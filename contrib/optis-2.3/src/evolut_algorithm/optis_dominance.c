/* Domination checking routines */
# include <math.h>

# include "../optis_global.h"

/* Routine for usual non-domination checking
   It will return the following values
   1 if a dominates b
   -1 if b dominates a
   0 if both a and b are non-dominated */

int check_dominance (individual *a, individual *b, int nobj)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;

    for (i=0; i<nobj; i++)
    {
	if(isnan(a->obj[i]))
	{
	    a->obj[i] = INF;
	}
 	if(isnan(b->obj[i]))
	{
	    b->obj[i] = INF;
	}


        if (a->obj[i] < b->obj[i])
        {
            flag1 = 1;

        }
        else
        {
            if (a->obj[i] > b->obj[i])
            {
                flag2 = 1;
            }
        }
    }
    if (flag1==1 && flag2==0)
    {
        return (1);
    }
    else
    {
        if (flag1==0 && flag2==1)
        {
            return (-1);
        }
        else
        {
            return (0);
        }
    }
}
