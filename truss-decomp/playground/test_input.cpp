//
// Created by yche on 10/10/19.
//

#include "util/file_system/file_util.h"
#include "util/log.h"
#include "util/program_options/popl.h"

int main(int argc, char *argv[]) {

    using namespace popl;
    OptionParser op("Allowed options");
    auto string_option = op.add<Value<std::string>>("f", "file-path", "the graph bin file path");

    op.parse(argc, argv);

    if (string_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);
        log_info("There are %zu triangles in the input graph.", 0);
    }
}