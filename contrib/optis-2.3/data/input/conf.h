/* Arquivo com as configurações necessárias para executar o programa.* 
 * Nao alterar a ordem dos dados.                                    */


/* Dados referentes ao algoritmo do NSGA-II. */
10				/* Número de gerações */
4				/* Tamanho da população (valor multiplo de 4) */
0.9				/* pcross_real, intervalo (0,1) */
0.5				/* pmut_real, intervalo (0,1) */
10				/* eta_c */
10				/* eta_m */
0.5				/* seed */


/* Informação do modelo calibrado */
/home/user/Programas/optis/data/input/struct.sh 			/* Script que gera a strutura do modelo */
inland            /* Binário do modelo (Dentro da estrutura gerada pelo script anterior) */
output/single_point-output.csv                    /* Arquivo de saída (Dentro da estrutura gerada pelo script anterior) */ 


/* Informações dos parâmetros (Dentro da estrutura gerada pelo script) */
4					/* Total de arquivos de parâmetros do modelo */
params/canopy				
params/fire				
params/soil				
params/vegetation			


/* Dados referentes aos sítios. */
2				/* Total de sítios */

km67				/* Nome do sítio */
26280			/* Total de dados comparados */
7		        	/* Número de colunas para o arquivo de dados observados */
13 			/* Número de colunas para o arquivo de output */
1 			/* Número de linhas a serem excluidas da leitura dos dados observados. As linhas excluidas podem ser um cabeçalho ou spin-up */
78841  			/* Número de linhas a serem excluidas da leitura dos outputs */
/home/user/Documentos/km67/obs_hourly.csv 		/* Nome do arquivo de dados observados */
3						/* Número de links simbolicos especificos de cada sítio */
/home/user/Documentos/km67/conf	data/offline/single_point/conf /* Diretório do conf do modelo para este sítio */
/home/user/Documentos/km67/input	data/offline/single_point/input /* Diretório do input do modelo */
/home/user/Documentos/km67/single_point_parameters data/offline/single_point/params/single_point_parameters			/* Nome do arquivo single_pointparameters do INLAND */

FNS				/* Nome do sítio */
26280 			/* Total de dados comparados */
7 			/* Número de colunas para o arquivo de dados observados */
13 		/* Número de colunas para o arquivo de output */
1 				/* Número de linhas a serem excluidas da leitura dos dados observados. As linhas excluidas podem ser um cabeçalho ou spin-up */
78841   			/* Número de linhas a serem excluidas da leitura dos outputs */
/home/user/Documentos/fns_km77_2/obs_hourly_fns.csv 		/* Nome dos arquivos de dados observados */
3
/home/user/Documentos/fns_km77_2/conf_fns data/offline/single_point/conf	/* Diretório do conf do modelo para este sítio */
/home/user/Documentos/fns_km77_2/input_fns data/offline/single_point/input	/* Diretório do input do modelo */
/home/user/Documentos/fns_km77_2/single_point_parameters_fns data/offline/single_point/params/single_point_parameters			/* Nome do arquivo single_pointparameters do INLAND */

/* Dados referentes aos níveis de calibração. Cada parâmetro deve ter seu limite máximo e mínimo. Deve-se informar para cada F.O. o seu nome, o sítio usado para calibrar essa variavel,      *
 * o tipo de medida de erro utilizada, a sua posição no arquivo de dados observados, e a sua posicao no arquivo de output */

2				/* Total de niveis hierárquicos */

1                               /* Modo de calibração (0: FO independente, 1:FO médias) */

Rn-H		/* Nome do nível de calibaração */
2								/* Total de parâmetros */
tau15 1 38 10 8 4000 5000                    			/* Nome, arquivo, linha, coluna de inicio, total de colunas, limite minimo e limite maximo do parâmetro */
ko15 1 40 10 8 0.1 0.5                    			/* Nome, arquivo, linha, coluna de inicio, total de colunas, limite minimo e limite maximo do parâmetro */

4				/* Total de funções objetivas */
Rn 1 0 4 5 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
H 1 0 5 8 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
Rn 2 0 4 5 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
H 2 0 5 8 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */

0                               /* Modo de calibaração (0: FO independente, 1:FO medias) */

LE-NEE2			/* Nome do nível de calibaração */
1								/* Total de parâmetros */
kc15 1 39 10 8 1.5e-05 2.5e-04                    			/* Nome, arquivo, linha, coluna de inicio, total de colunas, limite minimo e limite maximo do parâmetro */

4				/* Total de funções objetivas */
LE 1 0 6 7 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
NEE 1 0 7 10 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
LE 2 0 6 7 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
NEE 2 0 7 10 		/* Nome da variavel, sítio, medida de erro, pos observado e pos output */
