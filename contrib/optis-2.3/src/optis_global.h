# ifndef _GLOBAL_H_
# define _GLOBAL_H_

# include <stdio.h>

# define INF 1.0e14
# define EPS 1.0e-14

# define CALIB_MOD_IND 0
# define CALIB_MOD_MED 1

# define N_OF_METHODS 2
# define MAE 0
# define RMSE 1
# define MAE_STRING "MAE"
# define RMSE_STRING "RMSE"

# define ND_VALUE -9999.000

# define TEMPDIR "data/temp/"
# define B_SITIODIR "sitio"
# define AUXBASEDIR "aux"

# define INPUTDIR "data/input/"
# define OUTPUTDIR "data/output/"
# define OUTPARETODIR "data/output/paretoFront/"

# define CONF_FILE "conf.h"

/* Structs list */

typedef struct
{
    int rank;
    double *xreal;
    double *obj;
    double crowd_dist;
}individual;

typedef struct
{
    individual *ind;
}population;

typedef struct list
{
    int index;
    struct list *parent;
    struct list *child;
}list;

typedef struct
{
    double **obs_data;		/* Array com os dados observados */
    int nlines;		/* Numero de linhas de todos os arquivos (observado, flag e output) */
    int ncol_out;		/* Numero de colunas do arquivo de output */
    int exline_out;		/* Numero de linhas excluidas do inicio do arquivo de output (cabecalho ou spin-up) */
    int ncol_obs;		/* Numero de colunas do arquivo de dados observados */
    int exline_obs;		/* Numero de linhas excluidas do inicio do arquivo de observados (cabecalho ou spin-up) */
    char *name;			/* Nome do sitio */
    int nop;			/* Numero de ações (ex: links) para cada sítio */
    char ***op;			/* Operações a serem feitas */
}sitio;

typedef struct
{
    sitio *sit;
    int size;
}sit_array;

typedef struct
{
    int pos_obs; 		/* Posicao no arquivo de dados observados */
    int pos_out;		/* Posicao no arquivo de output do modelo */
    int pos_obj;		/* Posicao no array de funcoes objetivas de cada individuo */
}variavel;

typedef struct
{
    variavel *var;
    int size;
}var_array;

typedef struct
{
    int nfiles;			/* Numero de arquivos de parametros do modelo */
    char **filename;		/* Nome dos arquivos de parametros */
}base_param;

typedef struct
{
	int **cons;				/* First position the pft number, second the awood param, third aleaf and fourth aroot */
	int size;
}constraint_array;

typedef struct
{
    int cpars;			/* Total de parametros calibrados em cada nivel hierarquico */
    int *parl;		/* Array com a linha de cada parametro calibrado */
    int *parc;		/* Array com a coluna de cada parametro calibrado */
    int *pars;		/* Array com o tamanho de cada parametro calibrado */
    int *parf;		/* Array com o arquivo que contem cada parametro calibrado */
    char **parname;		/* Array com o tamanho de cada parametro calibrado */
    double *min_parvalue;	/* Array com o limite minimo para cada parametro calibrado */
    double *max_parvalue;	/* Array com o limite maximo para cada parametro calibrado */
    var_array **obj;		/* Array com 2 dim, onde 1 sao os sitios e 2 as medidas de erros */
    int nobj;			/* Total de funcoes objetivas do nivel hierarquico */
    int *obj_comp;			/* Array que indica por quantas variaveis cada objetivo e compartilhado (Versao multi-sitio com valor medio)*/
    int mode;			/* Modo de calibracao: 0 = calibracao de cada variavel indepedente, 1 = media de variaveis iguais */
    char *name;			/* Nome do nivel de calibracao */
    char **of_name;		/* Array com o nome de cada funcao objetiva, indicando a variavel e o metodo utilizado */
    int lvl;			/* Indica qual o nivel desta calibracao */
    constraint_array *constraint;	/* Contraint matriz to control the carbon allocation fraction (aroot + awood + aleaf = 1). */
}calibration_lvl;

typedef struct
{
    int ngen;			/* Total de geracoes do nsga */
    int popsize;		/* Tamanho da populacao */
    double pcross_real;
    double pmut_real;
    double eta_c;
    double eta_m;
    double seed;		/* Semente para a geracao de numeros randomicos */
}nsga_param;

typedef struct
{
    int ncal;			/* Total de calibracoes */
    double** min_value;	/* Array com o limite minimo para cada nivel de calibracao */
    double** max_value;	/* Array com o limite maximo para cada nivel de calibracao */
}paramsrange_opt;

typedef struct
{
    char *st_script;			/* nome do script que gera a estrutura necessária do modelo */
    char *bin;			    /* Nome do executavel do modelo */
    char *out_name;		/* Nome dos arquivos de output do modelo para cada escala temporal */
}modelo_st;

/* END - Structs list */


void run_inland(double *xreal, calibration_lvl *calib, base_param *bparams, int ind_id, int sitio_id, modelo_st*);

void write_inlandparams(char*, double*, base_param*, calibration_lvl*);

void store_params(individual *ind, calibration_lvl *calib, int popsize, base_param *bparams, sit_array *sitios);

double** rd_filewdouble(char*, int, int, int);

void evaluate_ind(int, individual*, calibration_lvl*, base_param*, sit_array*, modelo_st*);
void evaluate_pop(population*, int, calibration_lvl*, base_param*, sit_array*, modelo_st*);

double** rd_inlandout(int, int, int, int, int, modelo_st*);

void create_tempdir(int, sit_array*, modelo_st*);

individual* mindist_origin(population*, int, int);
double calcdist_origin(individual*, double*, int);


/* AUX functions */

char* rd_string(FILE*, int);
int rd_int(FILE*);
double rd_double(FILE*);

int** rd_filewint(char*, int, int, int);

nsga_param* init_nsgaparam(FILE*);

base_param* init_baseparam(FILE*);

sit_array* init_sitarray(FILE*);

int strcmpe(char*, char*);

void add_variavel(var_array*, variavel*);
calibration_lvl* init_calibrationlvl(FILE*, int, base_param*, sit_array*);
void rd_calibrationvars(int lvl, calibration_lvl* calib, FILE* file, sit_array* sitios);
void check_constraints(constraint_array *constraint, base_param *bparam);

char* rd_namewpath(FILE*, char*, int);

void set_paramsrange(int cal, int popsize, population *pop, int npars, paramsrange_opt *prange);
void init_paramsrange(paramsrange_opt *prange, int cal, int npars);

calibration_lvl** init_allcalibrationlvl(FILE *f_conf, base_param *bparam, sit_array *sitios, int ncal, calibration_lvl *allobj);

void add_objfromlvl(calibration_lvl*, calibration_lvl*, sit_array*);

modelo_st* init_modelo(FILE*);

void check_paramconstraint(char*, int, constraint_array* );

/* END - AUX functions */


/* Report functions */

FILE* open_allpopf(int, char*, int);
void report_pop(FILE*, population*, int, int, int);

FILE* open_bestindf();
void report_bestind(FILE*, individual*, calibration_lvl*);

FILE* open_calibparamsf();
void report_calibparams(calibration_lvl*, FILE*, individual*, char**);

void report_inlandparams(base_param*);

void report_ind(individual*, int, FILE*);

void report_paramsrange(paramsrange_opt *prange, calibration_lvl **allcaliblvl);

void report_allobj(calibration_lvl*, base_param*, sit_array*, modelo_st*);

void report_paretofront(base_param*, int, population*, calibration_lvl*);

/* END - Report functions */


/* Memory management */

population* alloc_pop(int, int, int);
void alloc_ind(individual*, int, int);

void deallocate2d(double**, int);

void dealloc_pop(population*, int);
void dealloc_ind(individual*);

nsga_param* alloc_nsgaparam();
void dealloc_nsgaparam(nsga_param*);

base_param* alloc_baseparam(int);
void dealloc_baseparam(base_param*);

sit_array* alloc_sitarray(int);
void alloc_sitio(sitio *sit);
void dealloc_sitarray(sit_array*);
void dealloc_sitio(sitio*);

void dealloc_calibrationlvl(calibration_lvl*, sit_array*);
constraint_array* alloc_constraint();
void dealloc_constraint(constraint_array *constraint);
void add_constraint(constraint_array *constraint);

calibration_lvl* alloc_allobj(sit_array *sitios);

paramsrange_opt* alloc_paramsrange(int ncal, calibration_lvl **allcaliblvl);
void dealloc_paramsrange(paramsrange_opt *prange);

void dealloc_allcalibrationlvl(int ncal, calibration_lvl **allcaliblvl, sit_array *sitios);

modelo_st* alloc_modelo();
void dealloc_modelo(modelo_st*);

/* END - Memory management */


/* Objective function methods */

double* calc_mae(double**, double**, int, var_array*);
double* calc_rmse(double**, double**, int, var_array*);

/* END - Objective function methods */


/* Evolutionary algorithm functions list. */

void randomize(double);
void warmup_random (double);
void advance_random (void);
double randomperc(void);
int rnd (int, int);
double rndreal (double, double);

void insert(list*, int);
list* del(list*);

int check_dominance(individual*, individual*, int);

void assign_rank_and_crowding_distance(population*, int, int);
void assign_crowding_distance_list(population*, list*, int, int);
void assign_crowding_distance_indices(population*, int, int, int);
void assign_crowding_distance(population*, int*, int**, int, int);

void quicksort_front_obj(population*, int, int[], int, int);
void quicksort_dist(population*, int*, int, int);

void selection(population*, population*, int, int, int, double, double*, double*, double);
individual* tournament(individual*, individual*, int);

void crossover(individual*, individual*, individual*, individual*, int, double, double*, double*, double);

void mutation_pop(population*, int, int, double, double*, double*, double);

void merge(population*, population*, population*, int, int, int);
void copy_ind(individual*, individual*, int, int);

void fill_nondominated_sort(population*, population*, int, int, int);
void crowding_fill(population*, population*, int, int, list*, int, int, int);

population* nsga(nsga_param*, calibration_lvl*, base_param*, sit_array*, modelo_st*);

void initialize_pop(population*,int, int, double*, double*);

/* END - Evolutionary algorithm funcitons. */

# endif
