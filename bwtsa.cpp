#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include "bwt_scan.hpp"

int main(int argc, char** argv) {
    std::string fname(argv[1]);
    std::vector<std::pair<uint64_t, uint64_t>> samples_first;
    std::vector<uint64_t> samples_last;
    bwt_scan_ssa(fname, samples_first, samples_last);
}
