# include <stdlib.h>
# include <string.h>
# include <sys/stat.h>

# include "optis_global.h"

void create_tempdir(int popsize, sit_array *sitios, modelo_st* modelo)
{
	char cmd[200];
	char dir[200], s_aux[200];
	int i, j, k;

	/* Deletando o diretorio TEMP antigo */
	sprintf(cmd, "rm -rf %s;", TEMPDIR);
	system(cmd);

	/* Criando um novo diretorio TEMP */
	mkdir(TEMPDIR, 0700);

	/* Criando diretorio para cada individuo, com todos os sitios */
	for(i=0;i<popsize;i++)
	{
		sprintf(dir, "%s%d/", TEMPDIR, i);
		mkdir(dir, 0700);

		for(j=0;j<sitios->size;j++)
		{
			sprintf(s_aux, "%s%s%d/", dir, B_SITIODIR, j);
			mkdir(s_aux, 0700);

			sprintf(cmd, "cd %s && %s;", s_aux, modelo->st_script);
			system(cmd);

			for(k=0;k<sitios->sit[j].nop;k++)
			{
				sprintf(cmd, "cd %s && ln -s %s %s;", s_aux, sitios->sit[j].op[k][0], sitios->sit[j].op[k][1]);
				system(cmd);

			}

		}
	}

	// Diretório auxiliar
	sprintf(dir, "%s%s/", TEMPDIR, AUXBASEDIR);
	mkdir(dir, 0700);

	sprintf(cmd, "cd %s && %s;", dir, modelo->st_script);
	system(cmd);

    return;
}
