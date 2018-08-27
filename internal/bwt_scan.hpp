#ifndef BWT_SCAN_H
#define BWT_SCAN_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "sdsl/wavelet_trees.hpp"

using std::cerr;
using std::cout;
using std::endl;

/* note: pfbwt uses 0x0 for end of string, 0x1 for end of word */
/* r-index uses 1 for end of string (we can change this, nbd?) */

template<typename S>
void get_F_and_runs(const S& str, std::vector<uint64_t>& F, sdsl::bit_vector& runstarts, sdsl::bit_vector& runends, uint64_t* term_pos=nullptr, uint8_t terminator=1) {
    /* NOTE: any char c s.t. c < terminator is treated as terminator (ie. if c==0 and terminator == 1, then c is treated as 1) */
    assert(runstarts.size() == str.size() && runends.size()==str.size() && F.size() == 256);
    F = vector<ulint>(256,0);
    uint8_t lastchar=str[0], nextchar=str[0], nextnextchar=str[1];
    uint64_t i;
    for (i = 0; i < str.size()-1; ++i) {
        nextnextchar = str[i+1];
        if (nextchar <= terminator) {
            F[terminator] += 1;
            if (term_pos) *term_pos = i;
        } else {
            F[nextchar] += 1;
        }
        if (nextchar != lastchar || i == 0) { // i is runstart or first char of bwt
            runstarts[i] = 1;
        }
        if (nextchar != nextnextchar) {
            runends[i] = 1;
        }
        lastchar = nextchar;
        nextchar = nextnextchar;
    }
    // deal with i==n-1 case (nextchar == str[n-2+1] == str[n-1], lastchar == str[n-2])
    if (nextchar <= terminator) F[terminator] += 1; else F[nextchar] += 1; 
    if (nextchar != lastchar || i == 0) runstarts[i] = 1;
    runends[i] = 1; // i==n-1 aka last character in bwt

    // make sure only one terminator exists
    if (F[terminator] > 1) {
        cerr << "error: bwt can only have one occurence of TERMINATOR! (found " << F[terminator] << " occurences" << endl;
        exit(1);
    }

    // properly finish populating F
    uint64_t running_sum = 1;
    uint64_t j;
    F[terminator] = 0;
    for (uint32_t i = terminator+1; i < F.size(); i++) {
        j = F[i];
        F[i] = running_sum;
        running_sum += j;
    }
}

void bwt_scan_ssa(std::string bwt_fname, std::vector<std::pair<uint64_t, uint64_t>>& samples_first, std::vector<uint64_t>& samples_last, std::vector<uint64_t>& F, uint64_t* term_pos=nullptr) {
    sdsl::wt_huff<> bwt;
    sdsl::construct(bwt, bwt_fname, 1);
    size_t n = bwt.size();
    sdsl::bit_vector run_starts(n,0);
    sdsl::bit_vector run_ends(n, 0);
    F.assign(256, 0);
    get_F_and_runs<sdsl::wt_huff<>>(bwt, F, run_starts, run_ends, term_pos);
    // rank-ify run_starts and run_ends
    sdsl::bit_vector::rank_1_type rs_rank(&run_starts);
    sdsl::bit_vector::rank_1_type re_rank(&run_ends);
    assert(run_starts.size() == run_ends.size());
    assert(sdsl::util::cnt_one_bits(run_starts) == sdsl::util::cnt_one_bits(run_ends));
    uint64_t nruns = sdsl::util::cnt_one_bits(run_starts);
    samples_first.resize(nruns);
    samples_last.assign(nruns, 0);
    uint64_t k = 0; // j is sa sample
    uint64_t j = n-1; // j is sa sample
    uint8_t c;
    // do LF mappings
    for (uint64_t i = 0; i < n; ++i) {
        // i is simply iteration no.
        // j is sa sample (starts at n)
        // k is bwt position (starts at 0 bc BWT[0]->T[n]
        if (run_starts[k]) {
            auto idx = rs_rank(k);
            samples_first[idx] = std::pair<uint64_t, uint64_t>(j ? j-1 : n-1, idx); 
        } 
        if (run_ends[k]) {
            auto idx = re_rank(k);
            samples_last[idx] = j ? j-1 : n-1;
        }
        // prepare for next iteration
        j--;
        k = F[bwt[k]] + bwt.rank(k, bwt[k]); // LF mapping, ie. next bwt position
    }
    /* DEBUG:
    std::ofstream saf("samples_first_bwtscan"); std::ofstream sal("samples_last_bwtscan");
    for (auto i: samples_first) saf << std::get<0>(i) << " " << std::get<1>(i) << endl;
    for (auto i : samples_last) sal << i << endl;
    saf.close(); sal.close();
    std::ofstream ffs("F");
    for (auto i = 0; i < F.size(); ++i) {
        ffs << i << " " << F[i] << endl;
    }
    */
}

#endif
