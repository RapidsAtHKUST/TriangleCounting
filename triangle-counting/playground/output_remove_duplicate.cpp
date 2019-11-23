//
// Created by yche on 11/23/19.
//

#include <cassert>
#include <malloc.h>

#include "ips4o/ips4o.hpp"

#include "util/pretty_print.h"
#include "util/popl.h"
#include "util/timer.h"
#include "util/util.h"
#include "util/log.h"
#include "primitives.h"
#include "yche_serialization.h"

using namespace std;
using namespace popl;
using namespace std::chrono;

#define PAGE_SIZE (4096)
#define IO_REQ_SIZE (PAGE_SIZE * 32)
#define IO_QUEUE_DEPTH (16)

inline std::string exec(const char *cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}


int main(int argc, char *argv[]) {
    OptionParser op("Allowed options");
    auto string_option = op.add<Value<std::string>>("f", "file-path", "the graph bin file path");
    auto output_string_option = op.add<Value<std::string>>("o", "output-path", "the ouptut bin file path");
    op.parse(argc, argv);

    using Edge = pair<int32_t, int32_t>;
    Timer global_timer;

    if (string_option->is_set() && output_string_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);

        auto file_name = string_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY | O_DIRECT, S_IRUSR | S_IWUSR);
        Edge *edge_lst = (Edge *) memalign(PAGE_SIZE, size + IO_REQ_SIZE);

#ifndef USE_LOG
        omp_set_num_threads(std::thread::hardware_concurrency());
#endif

        auto max_omp_threads = omp_get_max_threads();
        Timer io_timer;
        size_t read_size = 0;
#pragma omp parallel num_threads(IO_QUEUE_DEPTH)
        {
#pragma omp for schedule(dynamic, 1) reduction(+:read_size)
            for (size_t i = 0; i < size; i += IO_REQ_SIZE) {
                auto it_beg = i;
                auto *chars = reinterpret_cast<uint8_t *>(edge_lst);
                auto ret = pread(file_fd, chars + it_beg, IO_REQ_SIZE, it_beg);
                if (ret != IO_REQ_SIZE) {
                    log_error("Err, %zu, %zu, %zu, %d", i, it_beg, IO_REQ_SIZE, ret);
                } else {
                    read_size += ret;
                }
            }
#pragma omp single
            log_info("%zu, %zu", read_size, size);
        }
        log_info("IO Time: %.6lfs, DIO-QPS: %.6lf GB/s", io_timer.elapsed(), size / io_timer.elapsed() / pow(1024, 3));
        log_info("Load File Time: %.9lfs", global_timer.elapsed());

        // 1st: Remove Multi-Edges and Self-Loops.
        Timer sort_timer;
        int32_t max_node_id = 0;
#pragma omp parallel for reduction(max: max_node_id) schedule(dynamic, 32*1024)
        for (size_t i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
            }
            max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
        }
        log_info("Populate File Time: %.9lfs", global_timer.elapsed());

        ips4o::parallel::sort(edge_lst, edge_lst + num_edges, [](Edge l, Edge r) {
            if (l.first == r.first) {
                return l.second < r.second;
            }
            return l.first < r.first;
        });

#pragma omp parallel
        {
#pragma omp for
            for (size_t i = 0; i < num_edges - 1; i++) {
                auto l = edge_lst[i];
                auto r = edge_lst[i + 1];
                if (l.first == r.first) {
                    // = because of duplicates: self-loop-edge and multi-edge
                    assert(l.second <= r.second);
                } else {
                    assert(l.first < r.first);
                }
            }
#pragma omp single
            log_info("Verification Time: %.9lfs", sort_timer.elapsed());
        };
        log_info("Sort Time: %.9lfs", sort_timer.elapsed());
        auto num_vertices = static_cast<uint32_t >(max_node_id) + 1;
        log_info("Pre-Process Edge List Time: %.9lf s", global_timer.elapsed());

        auto graph_dir_path = output_string_option->value(0);
        exec(string("mkdir -p " + graph_dir_path).c_str());

        // Remove Duplicates Then.
        // Selection.
        Timer timer;
        vector<int32_t> histogram;
        auto edge_lst_buffer = (Edge *) memalign(PAGE_SIZE, size + IO_REQ_SIZE);
        auto *relative_off = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
#pragma omp parallel num_threads(max_omp_threads)
        {
            SelectNotFOMP(histogram, edge_lst_buffer, edge_lst, relative_off, num_edges, [edge_lst](uint32_t it) {
                return edge_lst[it].first == edge_lst[it].second || (it > 0 && edge_lst[it - 1] == edge_lst[it]);
            });
        }
        swap(edge_lst, edge_lst_buffer);
        num_edges = num_edges - relative_off[num_edges - 1];
        free(relative_off);
        log_info("New # of edges: %zu, Elapsed: %.9lfs", num_edges, timer.elapsed());
        log_debug("max_node_id: %d", max_node_id);

        FILE *pFile = fopen(string(graph_dir_path + "/" + "undir_edge_list.bin").c_str(), "w");
        YcheSerializer::write_array(pFile, edge_lst, num_edges);
        fflush(pFile);
        fclose(pFile);
    }
}