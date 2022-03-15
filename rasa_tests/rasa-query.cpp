// driver file to retrieve samples from an r-index using a combination of the tree
// data structure and phi.

#include <iostream>
#include <random>
#include <set>
#include <cstring>
#include <string>
#include "../internal/rle_string.hpp"
#include "../internal/r_index.hpp"
#include "../internal/utils.hpp"
#include "benchmark.hpp"
#include "rasa_common.hpp"

// #define DEBUG 1

using namespace ri;
using namespace std;

vector<ulint>& load_csa(std::string fname, ulint n, vector<ulint> &sa) {
    sa.clear();
    std::ifstream ifs(fname);
    uint64_t x = 0;
    uint64_t y = 0;
    while (ifs.read((char*) &x, 5) && ifs.read((char*) &y, 5)) {
        sa.push_back(y ? y-1 : n-1);
    }

    return sa;
}

int main(int argc, char** argv) {
  if(argc < 1){
    fprintf(stderr, "Usage: %s basename\n", argv[0]);
    return 1;
  }
  const string basename = argv[1];
  const string ri_filename = basename + ".ri";
  const string ssa_filename = basename + ".ssa";
  assert_file_exists(ri_filename.c_str());
  assert_file_exists(ssa_filename.c_str());

  r_index<> idx;
  idx.load_from_file(ri_filename.c_str());
  // string bwt = idx.get_bwt();
  // rle_string rle_bwt = rle_string(bwt);
  const auto text_length = idx.bwt.size();

  std::vector<ulint> ssa;
  load_csa(ssa_filename, text_length, ssa);
  std::cerr << "index, bwt, and rle loaded." << endl;
#ifdef DEBUG
  batch_query("rasa-index", basename, [&](uint64_t query_index) -> uint64_t { 
      my::Timer timer;
      timer.start();
      const auto comp_val = access_sa(idx, query_index);
      const double comp_time = timer.restart();
      const auto my_val = idx.query_csa(query_index, ssa);
      const double my_time = timer.restart();
      if(my_time < comp_time) {
        printf("SA[%lu]: my_time: %f comp_time : %f diff: %f\n", query_index, my_time, comp_time, my_time - comp_time);
      }
      CHECK_EQ(comp_val, my_val);
      return (my_val + 1) % text_length;
      }, true);
#else
  batch_query("rasa-index", basename, [&](uint64_t query_index) -> uint64_t { 
      return idx.query_csa(query_index, ssa);
      }, true);
#endif
  return 0;

  // ifstream answerfile(answer_filename);
  //
  // ifstream queryfile(query_filename);
  // while(true) {
  //   uint64_t answer_value;
  //   answerfile.read(reinterpret_cast<char*>(&answer_value), sizeof(decltype(answer_value)));
  //   uint64_t query_index;
  //   queryfile.read(reinterpret_cast<char*>(&query_index), sizeof(decltype(query_index)));
  //   if(!queryfile.good()) { break; }
  //    const uint64_t answer = idx.query_csa(query_index, ssa);
  //    DCHECK_EQ(answer, query_index);
  //    fwrite(&answer, sizeof(decltype(answer)), 1, stdout);
  // }
}
