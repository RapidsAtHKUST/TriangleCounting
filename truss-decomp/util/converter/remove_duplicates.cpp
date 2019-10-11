
//
// Created by yche on 10/11/19.
//

#include "util/program_options/popl.h"
#include "util/timer.h"
#include "util/log.h"
#include "util/sort/parasort_cmp.h"
#include "util/file_system/file_util.h"
#include "primitives.h"
#include "util/yche_serialization.h"

using namespace std;
using namespace popl;
using namespace std::chrono;

template<typename T, typename I>
T RemoveDuplicates(pair<T, T> *&edge_lst, I &num_edges) {
    using Edge = pair<T, T>;
    Timer timer;
    T max_node_id = 0;
    // Partition.
#pragma omp parallel
    {
#pragma omp for reduction(max: max_node_id)
        for (I i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
                max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
            }
        }
    }
    auto max_omp_threads = omp_get_max_threads();
    parasort(num_edges, edge_lst, [](const Edge &left, const Edge &right) {
        if (left.first == right.first) {
            return left.second < right.second;
        }
        return left.first < right.first;
    }, max_omp_threads);
    log_info("Finish Sort, %.9lfs", timer.elapsed());

    // Selection.
    auto *relative_off = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
    Edge *edge_lst2 = (Edge *) malloc(sizeof(Edge) * num_edges);
    auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);
#pragma omp parallel
    {
        FlagPrefixSumOMP(histogram, relative_off, num_edges, [edge_lst](uint32_t it) {
            return edge_lst[it].first == edge_lst[it].second || (it > 0 && edge_lst[it - 1] == edge_lst[it]);
        }, max_omp_threads);
#pragma omp for
        for (I i = 0u; i < num_edges; i++) {
            if (!(edge_lst[i].first == edge_lst[i].second || (i > 0 && edge_lst[i - 1] == edge_lst[i]))) {
                auto off = i - relative_off[i];
                edge_lst2[off] = edge_lst[i];
            }
        }
    }
    edge_lst = edge_lst2;
    num_edges = num_edges - relative_off[num_edges - 1];

    log_info("New # of edges: %zu, Elapsed: %.9lfs", num_edges, timer.elapsed());
    return max_node_id;
}

int main(int argc, char *argv[]) {
    OptionParser op("Allowed options");
    auto string_option = op.add<Value<std::string>>("f", "file-path", "the graph bin file path");
    auto output_option = op.add<Value<std::string>>("o", "output-path", "the output bin file path");
    op.parse(argc, argv);

    using Edge = pair<int32_t, int32_t>;
    Timer global_timer;
    if (string_option->is_set() && output_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);

        auto file_name = string_option->value(0);
        auto output_dir = output_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
        Edge *edge_lst = (Edge *) mmap(nullptr, size, PROT_READ | PROT_WRITE,
                                       MAP_PRIVATE | MAP_POPULATE, file_fd, 0);

        // Remove Multi-Edges and Self-Loops.
//        Edge *prev_edge_lst = edge_lst;
//        auto prev_num_edges = num_edges;
        auto max_node_id = RemoveDuplicates(edge_lst, num_edges);
        auto num_vertices = static_cast<uint32_t >(max_node_id) + 1;
        log_info("load edge list bin time: %.9lf s", global_timer.elapsed());

        string graphpath = output_dir + "/" + "undir_edge_list.bin";
        FILE *pFile = fopen(graphpath.c_str(), "w");
        YcheSerializer serializer;
        serializer.write_array(pFile, edge_lst, num_edges);
        fclose(pFile);
    }
}
