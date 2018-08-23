#include <cctype>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include <cctype>

#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread);

#include "internal/r_index.hpp"

#include "internal/utils.hpp"

using namespace ri;
using namespace std;
using namespace std::literals;

typedef r_index<> idx_t;

struct ri_opts_t {
    ri_opts_t() {};
    ulint max_hits = (ulint)-1;
    ulint max_range = (ulint)-1;
    uint z = 0;
};

struct sam_t {
    ~sam_t(){
        clear();
    }
    std::string query = "";
    uint flag = 4;
    std::vector<ulint> locs;
    std::string CIGAR = "";
    std::string seq = "";
    std::string qual = "";
    vector<std::string> tags;
    void clear() {
        query.clear();
        flag = 0;
        locs.clear();
        seq.clear();
        qual.clear();
        tags.clear();
    }
    friend std::ostream& operator <<(std::ostream&, const sam_t&);
};

std::ostream& operator <<(std::ostream& os, const sam_t& s) {
    if (s.locs.size()) {
        // for (const auto l: s.locs) {
        for (auto it = s.locs.begin(); it != s.locs.end(); ++it) {
            os << s.query << '\t';
            if (it == s.locs.begin()) os << s.flag << '\t';
            else os << (s.flag | 256) << '\t';
            os << "*" << '\t'; // ref: not implemented
            os << (*it)+1 << '\t';
            os << "*" << '\t'; // MAPQ: not implemented
            os << s.CIGAR << '\t';
            os << "*" << '\t'; // RNEXT
            os << "0" << '\t'; // PNEXT
            os << "0" << '\t'; // TLEN
            if (it == s.locs.begin()) os << s.seq << '\t' << s.qual;
            else os << "*" << '\t' << "*" ;
            if (s.tags.size()) {
                os << '\t';
                for (auto tag: s.tags) os << tag << " " ;
            }
            os << '\n';
        }
    } else {
        os << s.query << '\t' << s.flag << '\t' << "*" << '\t' << "*" << '\t' << "*" << '\t' << s.CIGAR << '\t' << "*" << '\t' << "0" << '\t' << "0" << '\t' << s.seq << '\t' << s.qual;
        if (s.tags.size()) {
            os << '\t'; 
            for (auto tag: s.tags) os << tag << " ";
        }
        os << '\n';
    }
    return os;
}


void help(){
    cout << "ri-align: align all occurrences of the input patterns." << '\n' << '\n';

    cout << "Usage: ri-align [options] <program> <index> <fastq>" << '\n';
    cout << "   <program>               locate/count/pw_locate         " << '\n';
    cout << "   <index>                 index file (with extension .ri)" << '\n';
    cout << "   <patterns>              file in pizza&chili format containing the patterns." << '\n';
    cout << "   --max-hits/-m <int>     maximum number of hits to return. [default: std::numeric_limits<size_t>::max()]\n";
    cout << "   --max-range/-r <int>  only return hits if a query has <= r occurences [default: std::numeric_limits<size_t>::max()]\n";
    cout << "   -z <int>                number of mismatches to allow\n";
    exit(0);
}


void parse_args(char** argv, int argc, int &ptr, size_t &niter, uint64_t &max_hits, uint64_t &max_range, uint32_t &z){

    assert(ptr<argc);

    string s(argv[ptr]);
    ptr++;

    if((s.compare("--max-hits") && s.compare("-m")) == 0) {
        if(ptr>=argc-1){
            cout << "Error: missing parameter after --max-hits/-m option." << '\n';
            help();
        }
        max_hits = std::strtoull(argv[ptr++], nullptr, 10);
    } else if ((s.compare("--max-range") && s.compare("-r")) ==0) {
        if(ptr>=argc-1){
            cout << "Error: missing parameter after --range-thres/-r option." << '\n';
            help();
        }
        max_range = std::strtoull(argv[ptr++], nullptr, 10);
    } else if(s.compare("-z") == 0) {
        if(ptr>=argc-1){
            cout << "Error: missing parameter after -z option." << '\n';
            help();
        }
        z = std::strtoul(argv[ptr++], nullptr, 10);
    } else{
        cout << "Error: unknown option " << s << '\n';
        help();
    }

}

/*
size_t count(idx_t& idx, kseq_t* seq, ri_opts_t opts, vector<sam_t>& sam) {
    sam.clear();
    ulint count = 0;
    if (opts.z == 0) {
        auto range = idx.exact_count(seq->seq.s, seq->seq.l);
        count = (range.second >= range.first) ? range.second - range.first + 1 : 0;
    } else {
        auto ranges = idx.inexact_count(seq->seq.s, seq->seq.l, opts);
        for (const auto r: ranges) {
            count += (r.second >= r.first) ? (r.second - r.first + 1) : 0;
        } 
    }
    cout << seq->name.s << "\t" << count << "\n";
    return count;
}

size_t locate(idx_t& idx, kseq_t* seq, ri_opts_t opts, vector<sam_t>& sams) {
    sam_t* sam = &(sams.at(0));
    sam->clear();
    if (opts.z == 0) {
        idx.exact_locate(std::string(seq->seq.s), opts, [sam, &seq](const std::string& P, const ri::range_t range, std::vector<ulint>& locs) {
                sam->clear();
                uint m = P.size();
                sam->query = std::string(seq->name.s);
                sam->CIGAR = std::to_string(m) + "M";
                sam->seq = std::string(seq->seq.s);
                sam->qual = std::string(seq->qual.s);
                sam->flag  = locs.size() ? 0 : 4;
                std::ostringstream tag_stream;
                // add tags here
                tag_stream << "NH:i:" << range.second-range.first+1;
                sam->tags.push_back(tag_stream.str());
                sam->locs.swap(locs);
                cout << *sam;
                });
    }
    return sam->locs.size();
}

size_t piecewise_count(idx_t& idx, kseq_t* seq, ri_opts_t opts, vector<sam_t>& sams) {
    (void) sams;
    idx.piecewise_count(std::string(seq->seq.s), opts, [&seq](const string& P, const ri::range_t range, const uint start, const uint end) {
            (void) P;
            (void) start;
            (void) end;
            cout << seq->name.s << "\t" << range.second - range.first + 1 << '\t' << end - start + 1 << endl;
            });
    return sams.size();
}

size_t piecewise_locate(idx_t& idx, kseq_t* seq, ri_opts_t opts, vector<sam_t>& sams) {
    sams.clear();
    idx.piecewise_locate(std::string(seq->seq.s), opts, [&seq, &sams](const string& P, const ri::range_t range, const uint start, const uint end, vector<ulint>& locs) {
            uint m = P.size();
            sam_t sam;
            sam.query = std::string(seq->name.s);
            sam.qual = std::string(seq->qual.s);
            sam.seq = std::string(seq->seq.s);
            if (end-start+1) {
                std::ostringstream CIGAR_stream;
                if (start) CIGAR_stream << start << "S";
                CIGAR_stream << end-start+1<< "M";
                if (m-end-1 > 0) CIGAR_stream << m-end-1 << "S";
                sam.CIGAR = CIGAR_stream.str();
            }
            if (!locs.size()) sam.flag = 4; // TODO: should we flag un-located alignments as 4? 
            else sam.flag = 0; 
            std::ostringstream tag_stream;
            // add tags here
            tag_stream << "NH:i:" << range.second-range.first+1;
            sam.tags.push_back(tag_stream.str());
            sam.locs.swap(locs);
            cout << sam;
            sams.push_back(std::move(sam));
            });
    // cerr << seq->name.s << "has " << sams.size() << " pieces." << endl;;
    return sams.size();
}
*/

size_t locate(idx_t& idx, kseq_t* seq, ri_opts_t opts, vector<sam_t>& sams) {
    std::string s(seq->seq.s);
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    vector<ulint> occs = idx.locate_all(s);
    return occs.size();
}

template<size_t (*F)(idx_t&, kseq_t*, ri_opts_t, vector<sam_t>&)>
void read_and_locate(std::string idx_pre, std::string patterns, size_t niter, ri_opts_t opts){
    std::string text;

    cerr << "loading fwd/rev r-index" << endl;
    // idx_t idx(string(idx_pre).append(".ri"), string(idx_pre).append(".rev.ri"));
    idx_t idx;
    std::ifstream ifs(string(idx_pre).append(".ri"));
    bool f;
    ifs.read((char*)&f, sizeof(f));
    idx.load(ifs);
    
    // LOAD IDX HERE

    cerr << "searching patterns ... " << '\n';

    /* read fastq here, search as we go*/

    std::fprintf(stderr, "Reading from fastq file\n");
    gzFile read_fp(gzopen(patterns.c_str(), "r"));
    if (read_fp == nullptr) {
        cerr << "invalid zip file" << '\n';
        exit(1);
    }
    kseq_t *seq(kseq_init(read_fp));
    vector<sam_t> sams;
    sam_t s;
    sams.push_back(s);

    size_t total_occs;
    size_t occs = 0;
    {
        vector<char> buf(1 << 16);
        std::sprintf(buf.data(), "Time: ");
        Timer t(buf.data());
        while (kseq_read(seq) >= 0) {
            occs = F(idx, seq, opts, sams);
            // cout << sams[0];
            total_occs += occs;
        }
    }
    kseq_destroy(seq);
    gzclose(read_fp);
}



int main(int argc, char** argv){

    if(argc < 3)
        help();

    int ptr = 1;
    size_t niter=1; 
    uint32_t z = 0;
    uint64_t max_hits=std::numeric_limits<uint64_t>::max();
    uint64_t max_range=std::numeric_limits<uint64_t>::max();

    while(ptr<argc-3)
        parse_args(argv, argc, ptr, niter, max_hits, max_range, z);

    string program(argv[ptr]);
    string idx_pre(argv[ptr+1]);
    string patt_file(argv[ptr+2]);

    // std::ifstream in(pre);

    // bool fast;

    //fast or small index?

    ri_opts_t opts;
    opts.max_hits = max_hits;
    opts.max_range = max_range;
    opts.z = z;
    cerr << "Loading r-index (parameters: " << max_hits << ", " << max_range << ", " << z << ")\n";
    if (program == "locate") {
        read_and_locate<locate>(idx_pre, patt_file, niter, opts);
    } /* else if (program == "count") {
        read_and_locate<count>(idx_pre, patt_file, niter, opts);
    } else if (program == "pw_locate") {
        read_and_locate<piecewise_locate>(idx_pre, patt_file, niter, opts);
    } else if (program == "pw_count") {
        read_and_locate<piecewise_count>(idx_pre, patt_file, niter, opts);
    }
    */
}
