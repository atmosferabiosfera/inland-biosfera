# include "../optis_global.h"

FILE* open_allpopf(int clvl, char *cnome, int nobj)
{
    FILE *file;
    char name[500];

    sprintf(name, "%s%2.2d.%s.all_pop.out", OUTPUTDIR, clvl, cnome);
    file = fopen(name, "w");

    fprintf(file, "# This file contains the data of all generations\n");
    fprintf(file, "# Each individual contains %d objectives, rank, crowding distance\n", nobj);

    return file;
}
