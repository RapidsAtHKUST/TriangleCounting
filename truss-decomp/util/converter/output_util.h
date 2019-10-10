#pragma once

#include <string>
#include <vector>

using namespace std;

std::string exec(const char *cmd);

void WriteToOutputFiles(string &deg_output_file, string &adj_output_file, vector<int> &degrees,
                        vector<int32_t> &dst_vertices);

