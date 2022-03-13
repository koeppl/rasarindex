#include "../internal/rle_string.hpp"
#include "../internal/r_index.hpp"
#include "../internal/utils.hpp"

ulint access_sa(ri::r_index<>& idx, uint query_index) {
  auto& rle_bwt = idx.bwt;
  const ulint run = rle_bwt.run_of_position(query_index); // run number
  const ulint j = rle_bwt.run_range(run).second; // run ends at this value
  ulint phi_val = idx.samples_last[run]+1; // end sample at specified run

  for (size_t iter = 0; iter < (j - query_index); iter++) {
    phi_val = idx.Phi(phi_val); // sa-1 until we find n (j-i times)
  }
  return phi_val;
}
