# include <string.h>
# include <ctype.h>

# include "../optis_global.h"

int strcmpe(char *s1, char *s2)
{
	int i, size;
	char c1[100], c2[100];

	strcpy(c1, s1);
	strcpy(c2, s2);

	size = strlen(c1);
	for(i=0;i<size;i++)
	{
		c1[i] = toupper(c1[i]);
	}

	size = strlen(c2);
	for(i=0;i<size;i++)
	{
		c2[i] = toupper(c2[i]);
	}

	return strcmp(c1, c2);
}
