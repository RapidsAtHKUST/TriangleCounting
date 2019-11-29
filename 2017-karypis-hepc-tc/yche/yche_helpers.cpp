//
// Created by yche on 5/20/19.
//

extern "C"
{
#include <GKlib.h>
#include "../struct.h"
#include "../defs.h"
}

#include "pre_processing.h"
//#define LOCKS // Problematic, do not know why yet.
#ifdef LOCKS

#include <omp.h>

#endif

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

#include "yche_serialization.h"
#include "log.h"
#include "timer.h"


using namespace std;
using namespace std::chrono;

using edge_lst_type = vector<pair<int32_t, int32_t>>;

edge_lst_type ReadEdgeList(string dir) {
    log_info("Path: %s", dir.c_str());
    edge_lst_type edges;
    string graphpath = dir + "/" + "undir_edge_list.bin";
    FILE *pFile = fopen(graphpath.c_str(), "r");
    YcheSerializer serializer;
    serializer.read_array(pFile, edges);
    fclose(pFile);
    log_info("edge#: %zu", edges.size());
    return edges;
}

void ConvertEdgeListToCSRWrapper(gk_graph_t *&graph, vector<pair<int32_t, int32_t>> &edge_lst) {
    // Assume no duplicates.
    int32_t max_node_id = 0;
#pragma omp parallel for reduction(max: max_node_id) schedule(dynamic, 32*1024)
    for (size_t i = 0u; i < edge_lst.size(); i++) {
        max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
    }

    // at most 2G vertices
    graph = gk_graph_Create();
    graph->nvtxs = max_node_id + 1;
    graph->xadj = gk_zsmalloc(graph->nvtxs + 1, 0, "xadj");

    uint32_t *deg_lst;
    size_t *off;
    ConvertEdgeListToCSR(edge_lst.size(), &edge_lst.front(), graph->nvtxs, deg_lst, off,
                         graph->adjncy, omp_get_max_threads());
    free(deg_lst);

    memcpy(graph->xadj, off, sizeof(ssize_t) * (graph->nvtxs + 1));
    free(off);
}

#ifdef __cplusplus
extern "C" {
#endif

gk_graph_t *gk_graph_Read_bin_edge_lst(char *filename) {
    Timer timer;

    auto edge_lst = ReadEdgeList(filename);
    log_info("load edge list bin time: %.3lf s", timer.elapsed());
    gk_graph_t *graph = gk_graph_Create();
    ConvertEdgeListToCSRWrapper(graph, edge_lst);
    return graph;
}

#ifdef __cplusplus
}
#endif
