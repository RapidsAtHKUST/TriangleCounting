//
// Created by yche on 10/10/19.
//

#include "util/file_system/file_util.h"
#include "util/log.h"
#include "util/program_options/popl.h"

int main(int argc, char *argv[]) {

    using namespace popl;
    using namespace std;
    OptionParser op("Allowed options");
    auto string_option = op.add<Value<std::string>>("f", "file-path", "the graph bin file path");

    op.parse(argc, argv);
    using Edge = pair<uint32_t, uint32_t>;

    if (string_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);

        auto file_name = string_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
        Edge *edge_list =
                (Edge *) mmap(0, size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, file_fd, 0);
        log_info("There are %zu triangles in the input graph.", 0);
    }
}