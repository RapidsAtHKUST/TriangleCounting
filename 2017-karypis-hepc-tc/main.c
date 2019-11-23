/*!
\file
\brief The entry point of the triangle counting code
\date Started 5/10/2017
\author George
\version\verbatim $Id: cmdline.c 20946 2017-05-10 23:12:48Z karypis $ \endverbatim
*/

#include "tc.h"
#include <yche/log.h>

/*************************************************************************
* The entry point 
**************************************************************************/
int main(int argc, char *argv[]) {
    int64_t ntriangles = 0;
    params_t *params;
    vault_t *vault;

    params = getcmdline_params(argc, argv);

    vault = params->vault = loadData(params);
    vault->nprobes = 0;

#if defined(_OPENMP)
    omp_set_num_threads(params->nthreads);
    omp_set_nested(0);
    omp_set_dynamic(0);
#endif

    log_info("-----------------");
    log_info("  infile: %s", params->infile);
    log_info("  #nvtxs: %d", vault->graph->nvtxs);
    log_info(" #nedges: %zd", vault->graph->xadj[vault->graph->nvtxs]);
#if defined(_OPENMP)
    log_info("nthreads: %d", omp_get_max_threads());
#else
    params->nthreads = 1;
    printf("nthreads: 1 (was not compiled with openmp support)\n");
#endif
    log_info("");


    gk_startwctimer(vault->timer_global);
    ntriangles = ptc_MapJIK(params, vault);
    gk_stopwctimer(vault->timer_global);

    log_info("Results...");
    log_info("  #triangles: %12"PRId64"; #probes: %12"PRIu64"; rate: %10.2lf MP/sec",
           ntriangles, vault->nprobes,
           ((double) vault->nprobes) / ((double) 1e6 * gk_getwctimer(vault->timer_tc)));

    log_info("Timings...");
    log_info("     preprocessing: %9.3lfs", gk_getwctimer(vault->timer_pp));
    log_info(" triangle counting: %9.3lfs", gk_getwctimer(vault->timer_tc));
    log_info("    total (/x i/o): %9.3lfs", gk_getwctimer(vault->timer_global));
    log_info("-----------------");
}

