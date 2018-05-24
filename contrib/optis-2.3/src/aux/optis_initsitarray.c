# include <stdlib.h>
# include <string.h>
# include <unistd.h>

# include "../optis_global.h"

sit_array* init_sitarray(FILE *file)
{
    sit_array *sitios;
    int nsitio, i, j;
    char *f_name;

    nsitio = rd_int(file);

    if(nsitio < 1){
    	printf("\nERRO: O total de sitios precisa ser maior que zero.\n");
    	exit(EXIT_FAILURE);
    }

    sitios = alloc_sitarray(nsitio);

    for(i=0;i<nsitio;i++)
    {
    	/* Lendo nome do sitio */
    	sitios->sit[i].name = rd_string(file, 100);

    	/* Lendo total de linhas dos arquivos */
    	sitios->sit[i].nlines = rd_int(file);
    	if(sitios->sit[i].nlines < 0){
    		printf("\nERRO: O total de dados no sitio %d precisa ser positivo.\n", i+1);
    		exit(EXIT_FAILURE);
    	}

    	/* Lendo o numero de colunas do arquivo de dados observados */
    	sitios->sit[i].ncol_obs = rd_int(file);
    	if(sitios->sit[i].ncol_obs < 0){
    		printf("\nERRO: O total de colunas do arquivo de dados observados no sitio %d precisa ser positivo.\n", i+1);
    		exit(EXIT_FAILURE);
    	}

    	/* Lendo o numero de coluna do arquivo de output */
   		sitios->sit[i].ncol_out = rd_int(file);
    	if(sitios->sit[i].ncol_out < 0){
    		printf("\nERRO: O total de colunas do arquivo de output no sitio %d precisa ser positivo.\n", i+1);
    		exit(EXIT_FAILURE);
    	}

    	/* Lendo o numero de linhas a serem excluidas do arquivo de dados observados */
   		sitios->sit[i].exline_obs = rd_int(file);
   		if(sitios->sit[i].exline_obs < 0){
   			printf("\nERRO: O total de linhas excluidas do arquivo de dados observados no sitio %d precisa ser positivo.\n", i+1);
   			exit(EXIT_FAILURE);
   		}

    	/* Lendo o numero de linhas a serem excluidas do output */
   		sitios->sit[i].exline_out = rd_int(file);
   		if(sitios->sit[i].exline_out < 0){
   			printf("\nERRO: O total de linhas excluidas do arquivo de output no sitio %d precisa ser positivo.\n", i+1);
   			exit(EXIT_FAILURE);
   		}

    	/* Lendo os dados observados */
   		f_name = rd_string(file, 300);
		if(access(f_name, R_OK) != 0){
			printf("\nERRO: O arquivo %s nao existe ou nao pode ser lido.\n", f_name);
    		exit(EXIT_FAILURE);
    	}
    	sitios->sit[i].obs_data = rd_filewdouble(f_name, sitios->sit[i].nlines, sitios->sit[i].ncol_obs, sitios->sit[i].exline_obs);
   		free(f_name);

    	sitios->sit[i].nop = rd_int(file);

    	sitios->sit[i].op = (char***)malloc(sitios->sit[i].nop*sizeof(char**));

    	/* Lendo as operações que serão feitas */
    	for(j=0;j<sitios->sit[i].nop;j++)
    	{
    		sitios->sit[i].op[j] = (char**)malloc(2*sizeof(char*));

    		sitios->sit[i].op[j][0] = rd_string(file, 300);
    		sitios->sit[i].op[j][1] = rd_string(file, 300);

    		if(access(sitios->sit[i].op[j][0], R_OK) != 0){
    		    printf("\nERRO: O arquivo ou diretorio %s nao existe.\n", sitios->sit[i].op[j][0]);
    		    exit(EXIT_FAILURE);
    		}
    	}

    }


    return sitios;
}
