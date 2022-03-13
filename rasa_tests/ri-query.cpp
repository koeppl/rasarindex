// driver file to query an r-index using only phi.

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

using namespace ri;
using namespace std;


int main(int argc, char** argv) {
  if(argc < 1){
    fprintf(stderr, "Usage: %s basename\n", argv[0]);
    return 1;
  }
  const string basename = argv[1];
  const string ri_filename = basename + ".ri";
  assert_file_exists(ri_filename.c_str());

  r_index<> idx;
  idx.load_from_file(ri_filename.c_str());

  // string bwt = idx.get_bwt();
  // rle_string rle_bwt = rle_string(bwt);
  const auto text_length = idx.bwt.size();


  std::cerr << "index, bwt, and rle loaded." << endl;
  batch_query("r-index", basename, [&](uint64_t query_index) -> uint64_t { 
  return (access_sa(idx, query_index) + 1) % text_length;
      }
      , false);
}
