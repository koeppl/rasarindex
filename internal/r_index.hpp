/*
 * r_index.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Main class implementing the r-index
 *
 */

#ifndef R_INDEX_H_
#define R_INDEX_H_

#include <definitions.hpp>
#include <rle_string.hpp>
#include <sparse_sd_vector.hpp>
#include <permutation.hpp>

using namespace sdsl;

namespace ri{

class r_index{

public:

	using rle_string_t	= rle_string_sd;	//run-length encoded string

	r_index(){}

	/*
	 * Build index
	 */
	r_index(string &input){

		assert(not contains_reserved_chars(input));

		//build run-length encoded BWT
		{

			string bwt_s = build_bwt(input);
			bwt = rle_string_t(bwt_s);

			//build F column
			F = vector<ulint>(256,0);
			for(uchar c : bwt_s)
				F[c]++;

			for(ulint i=255;i>0;--i)
				F[i] = F[i-1];

			F[0] = 0;

			for(ulint i=1;i<256;++i)
				F[i] += F[i-1];

			for(ulint i=0;i<bwt_s.size();++i)
				if(bwt_s[i]==TERMINATOR)
					terminator_position = i;

			assert(input.size()+1 == bwt.size());

		}

		auto r = bwt.number_of_runs();

		//mark positions on bitvectors U,D
		//U, D have size bwt.size(): we include text terminator

		/*
		 * this vector stores the text position associated to the last BWT position
		 * of each BWT run, except the last.
		 */
		vector<ulint> SA_down(r-1);

		{

			vector<bool> U_vec(bwt.size(),false);
			vector<bool> D_vec(bwt.size(),false);

			//current BWT position (in first cycle will become terminator_position)
			ulint k = bwt.size();

			//current text position corresponding to bwt position k
			ulint j = bwt.size()-1;

			//repeat until we are back to the terminator
			while(k != terminator_position){

				if(k==bwt.size()) k = terminator_position;

				//if not first or last char in the BWT
				if(k > 0 and k < bwt.size()-1){

					assert(not U_vec[j]);
					assert(not D_vec[j]);

					//in this case, k is the first position of its run (Up position)
					if(bwt[k] != bwt[k-1]){

						U_vec[j] = true;

					}

					//in this case, k is the last position of its run (Down position)
					if(bwt[k] != bwt[k+1]){

						D_vec[j] = true;

						assert(bwt.run_of_position(k)<r-1);
						SA_down[bwt.run_of_position(k)] = j;

					}

				}

				k = LF(k);
				assert(k == terminator_position or j>0);
				j--;

			}

			//build gap-encoded bitvectors
			U = sparse_sd_vector(U_vec);
			D = sparse_sd_vector(D_vec);

			cout << D.rank(D.size()) << " " << r << endl;
			cout << U.rank(U.size()) << " " << r << endl;

			assert(D.rank(D.size())==r-1);
			assert(U.rank(U.size())==r-1);

		}


		//now build the invertible permutations UD, DR

		{

			vector<ulint> UD_vec(r-1,r);
			vector<ulint> DR_vec(r-1,r);

			//current BWT position (in first cycle will become terminator_position)
			ulint k = bwt.size();

			//current text position corresponding to bwt position k
			ulint j = bwt.size()-1;

			//repeat until we are back to the terminator
			while(k != terminator_position){

				if(k==bwt.size()) k = terminator_position;

				//if not first or last char in the BWT
				if(k > 0 and k < bwt.size()-1){

					//in this case, k is the first position of its run (Up position)
					if(bwt[k] != bwt[k-1]){

						assert(D[SA_down[ bwt.run_of_position(k-1)]]);
						assert(U[j]);

						ulint down = D.rank( SA_down[ bwt.run_of_position(k-1)] );
						ulint up = U.rank( j );

						//k must be in the next run of k-1
						assert(bwt.run_of_position(k) == bwt.run_of_position(k-1)+1);

						assert(up<r-1);
						assert(down<r-1);
						UD_vec[up] = down;

					}

					//in this case, k is the last position of its run (Down position)
					if(bwt[k] != bwt[k+1]){

						assert(D[j]);
						assert(bwt.run_of_position(k)<r-1);

						DR_vec[ D.rank(j) ] = bwt.run_of_position(k);

					}

				}

				k = LF(k);
				assert(k == terminator_position or j>0);
				j--;

			}

			assert(not_contains(UD_vec,r));
			assert(not_contains(DR_vec,r));

			UD = permutation(UD_vec);
			DR = permutation(DR_vec);

		}

	}

	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//inclusive range
		return {0,bwt_size()-1};

	}

	uchar operator[](ulint i ){
		return bwt[i];
	}

	/*
	 * \param r inclusive range of a string w
	 * \param c character
	 * \return inclusive range of cw
	 */
	range_t LF(range_t rn, uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		//number of c before the interval
		ulint c_before = bwt.rank(rn.first,c);

		//number of c inside the interval rn
		ulint c_inside = bwt.rank(rn.second+1,c) - c_before;

		//if there are no c in the interval, return empty range
		if(c_inside==0) return {1,0};

		ulint l = F[c] + c_before;

		return {l,l+c_inside-1};

	}

	//backward navigation of the BWT
	ulint LF(ulint  i){

		auto c = bwt[i];
		return F[c] + bwt.rank(i,c);

	}

	//forward navigation of the BWT
	ulint FL(ulint  i){

		//i-th character in first BWT column
		auto c = F_at(i);

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	//forward navigation of the BWT, where for efficiency we give c=F[i] as input
	ulint FL(ulint  i, uchar c){

		//i-th character in first BWT column
		assert(c == F_at(i));

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	/*
	 * access column F at position i
	 */
	uchar F_at(ulint i){

		ulint c = (upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
		assert(c<256);
		assert(i>=F[c]);

		return uchar(c);

	}

	/*
	 * Return BWT range of character c
	 */
	range_t get_char_range(uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		ulint l = F[c];
		ulint r = bwt_size()-1;

		if(c<255)
			r = F[c+1]-1;

		return {l,r};

	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){

		auto range = full_range();
		ulint m = P.size();

		for(ulint i=0;i<m and range.second>=range.first;++i)
			range = LF(range,P[m-i-1]);

		return range;

	}

	/*
	 * iterator locate(string &P){
	 *
	 * 	return iterator to iterate over all occurrences without storing them
	 * 	in memory
	 *
	 * }
	 */

	/*
	 * get number of runs in the BWT (terminator character included)
	 */
	ulint number_of_runs(){
		return bwt.number_of_runs();
	}

	/*
	 * get terminator (0x1) position in the BWT
	 */
	ulint get_terminator_position(){
		return terminator_position;
	}

	/*
	 * get BWT in string format
	 */
	string get_bwt(){
		return bwt.toString();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint w_bytes = 0;

		assert(F.size()>0);
		assert(bwt.size()>0);

		out.write((char*)&terminator_position,sizeof(terminator_position));
		out.write((char*)F.data(),256*sizeof(ulint));

		w_bytes += sizeof(terminator_position) + 256*sizeof(ulint);

		w_bytes += bwt.serialize(out);

		w_bytes += U.serialize(out);
		w_bytes += D.serialize(out);
		w_bytes += UD.serialize(out);
		w_bytes += DR.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&terminator_position,sizeof(terminator_position));

		F = vector<ulint>(256);
		in.read((char*)F.data(),256*sizeof(ulint));

		bwt.load(in);

		U.load(in);
		D.load(in);
		UD.load(in);
		DR.load(in);

	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".ri" will be automatically added
	 */
	void save_to_file(string path_prefix){

		string path = string(path_prefix).append(".ri");

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * load the structure from the path specified.
	 * \param path: full file name
	 */
	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}

	ulint text_size(){
		return bwt.size()-1;
	}

	ulint bwt_size(){
		return bwt.size();
	}

	uchar get_terminator(){
		return TERMINATOR;
	}

	string get_extension(){
		return string(".ri");
	}

	ulint print_space(){

		cout << "Number of runs = " << bwt.number_of_runs() << endl<<endl;

		ulint tot_bytes = bwt.print_space();

		cout << "\nTOT BWT space: " << tot_bytes << " Bytes" <<endl<<endl;

		return tot_bytes;

	}

private:

	bool not_contains(vector<ulint> &V, ulint x){

		ulint r=0;

		for(auto y:V){

			if(y==x){

				cout << "failed at run " << r << " / " << V.size() << endl;
				return false;

			}

			r++;

		}

		return true;

	}

	uint8_t bitsize(uint64_t x){

		if(x==0) return 1;
		return 64 - __builtin_clzll(x);

	}

	/*
	 * check if range rn on column F contains a
	 * single character
	 */
	bool uniq_char(range_t rn){

		for(ulint i=0;i<256;++i){

			ulint l = F[i];
			ulint r = (i==255?bwt_size():F[i+1]);

			if(rn.first >= l and rn.second < r ) return true;

		}

		return false;

	}

	/*
	 * builds BWT of input string using SE_SAIS algorithm
	 * uses 0x1 character as terminator
	 *
	 */
	static string build_bwt(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = SE_SAIS;
	    construct_sa<8>(cc);

	    //now build BWT from SA
	    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cc));

	    {

	        for (ulint i=0; i<sa.size(); i++){
	            auto x = sa[i];

	            assert(x<=text.size());

	            if ( x > 0 )
	            	bwt_s.push_back((uchar)text[x-1]);
	            else
	            	bwt_s.push_back(TERMINATOR);

	        }

	    }

	    sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
	    sdsl::remove(cache_file_name(conf::KEY_SA, cc));

	    return bwt_s;

	}

	static bool contains_reserved_chars(string &s){

		for(auto c : s)
			if(c == 0 or c == 1)
				return true;

		return false;

	}

	static const uchar TERMINATOR = 1;


	/*
	 * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
	 */

	//F column of the BWT (vector of 256 elements)
	vector<ulint> F;
	//L column of the BWT, run-length compressed
	rle_string_t bwt;
	ulint terminator_position = 0;

	/*
	 * Invertible permutations (i.e. efficient map and inverse):
	 *
	 * UD maps text positions associated to first (up) position of each BWT run to
	 * the last (down) position of the previous run. Size of the permutation is r-1 because
	 * the last run does not have runs after it.
	 *
	 * DR  maps text positions associated to Last (down) position of each BWT run to
	 * the corresponding BWT run (i.e. last position of that run)
	 *
	 * total: 2r log r * (1+epsilon) bits
	 *
	 */

	permutation UD; //Up samples to Down samples.
	permutation DR; //Down samples to BWT runs (last position of each run). Last run is excluded!

	/*
	 * gap-encoded bitvectors marking with a bit set sampled positions on the text
	 *
	 * U : marks text positions that are the first in their BWT run, except for the first run.
	 * D : marks text positions that are the last in their BWT run, except for the last run.
	 *
	 * r-1 bits set in each bitvector. Overall size = 2r log(n/r) bits
	 *
	 */

	sparse_sd_vector U;
	sparse_sd_vector D;

	/*
	 * overall: UD, DR, U, and D take r log n bits (plus low-order terms)
	 */

};

}

#endif /* R_INDEX_H_ */
