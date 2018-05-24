# include "../optis_global.h"

population* nsga(nsga_param *nparam, calibration_lvl *calib, base_param *bparam, sit_array *sitios, modelo_st *modelo)
{
    int i;
    FILE *apfile;
    int nobj, nreal, ngen, popsize;
    double *min_realvar, *max_realvar, pcross_real, pmut_real, eta_c, eta_m, seed;
    population *parent_pop, *child_pop, *mixed_pop;

    nobj = calib->nobj;
    nreal = calib->cpars;
    min_realvar = calib->min_parvalue;
    max_realvar = calib->max_parvalue;
    ngen = nparam->ngen;
    popsize = nparam->popsize;
    pcross_real = nparam->pcross_real;
    pmut_real = nparam->pmut_real;
    eta_c = nparam->eta_c;
    eta_m = nparam->eta_m;
    seed = nparam->seed;

    parent_pop = alloc_pop(popsize, nreal, nobj);
    child_pop = alloc_pop(popsize, nreal, nobj);
    mixed_pop = alloc_pop(2*popsize, nreal, nobj);

    randomize(seed);

    apfile = open_allpopf(calib->lvl, calib->name, nobj);

    initialize_pop (parent_pop, popsize, nreal, min_realvar, max_realvar);
    evaluate_pop (parent_pop, popsize, calib, bparam, sitios, modelo);
    assign_rank_and_crowding_distance (parent_pop, popsize, nobj);
    report_pop(apfile, parent_pop, 1, popsize, nobj);
    printf("\n gen = 1, done!\n\n");

    for (i=2; i<=ngen; i++)
    {
        selection (parent_pop, child_pop, popsize, nobj, nreal, pcross_real, min_realvar, max_realvar, eta_c);
        mutation_pop (child_pop, popsize, nreal, pmut_real, min_realvar, max_realvar, eta_m);
        evaluate_pop(child_pop, popsize, calib, bparam, sitios, modelo);
        merge (parent_pop, child_pop, mixed_pop, popsize, nreal, nobj);
        fill_nondominated_sort (mixed_pop, parent_pop, popsize, nobj, nreal);
        report_pop(apfile, parent_pop, i, popsize, nobj);
        printf("\n gen = %d, done!\n\n",i);
    }


    dealloc_pop(child_pop, popsize);
    dealloc_pop(mixed_pop, 2*popsize);

    fclose(apfile);
    
    printf("\n NSGA-II routine successfully exited \n");

    return parent_pop;
}
