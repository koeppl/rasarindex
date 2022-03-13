#ifndef R_INDEX_R_A_
#define R_INDEX_R_A_

#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <definitions.hpp>
#include <r_index.hpp>
#include <ri_rasa_tree.hpp>
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

using namespace sdsl;

// temp notes for me:
// - pred_to_run is basically sa_map so that means it can be removed
// - we need to do < for checking the costs because of the fact that we are counting them along the way.

namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>

//! Builds the graph of the SA. Finds cycles within the SA. Then builds and stores the trees built on them.
class rads {
public:
  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    build_rads_phi(unsorted_ssa, ssa, esa, pred);
    cout << "Searching paths ..." << endl;
    find_cycles();
    cout << "\nDebug info: " << endl;
    print_debug_info();
  }

  rads(const rads &other_rads) {
    this->sa_map = other_rads.sa_map;
    this->phi_x_inv = other_rads.phi_x_inv;
    this->trees = other_rads.trees;
    this->tree_pointers = other_rads.tree_pointers;
    this->bounds = other_rads.bounds;
    this->sa_graph = other_rads.sa_graph;
    this->trees_bv = other_rads.trees_bv;
  }

  rads(rads &&other_rads)
  : sa_map(move(other_rads.sa_map))
  , phi_x_inv(move(other_rads.phi_x_inv))
  , trees(move(other_rads.trees))
  , tree_pointers(move(other_rads.tree_pointers))
  , bounds(move(other_rads.bounds))
  , sa_graph(move(other_rads.sa_graph))
  , trees_bv(move(other_rads.trees_bv))
  {}

  rads& operator=(const rads &other_rads) {
    return *this = rads(other_rads);
  }

  rads& operator=(rads &&other_rads) {
    swap(sa_map, other_rads.sa_map);
    swap(phi_x_inv, other_rads.phi_x_inv);
    swap(trees, other_rads.trees);
    swap(tree_pointers, other_rads.tree_pointers);
    swap(bounds, other_rads.bounds);
    swap(sa_graph, other_rads.sa_graph);
    swap(trees_bv, other_rads.trees_bv);
    return *this;
  }

  //! Build the csa/rads to be used with Phi^-1.
  /*!
    \param unsorted_ssa Non-sorted start samples. -> sorted in SA order
    \param ssa Sorted start samples.
    \param esa Non-sorted end samples.
    \param pred Predecessor bit vector. Used to calculate successors and predecessors.
  */
  void build_rads_phi_inv(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa;
    sa_map.reserve(esa.size()); //! TODO: not needed!
    bounds.resize(esa.size());
    sa_graph.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[esa[i]] = i; // sa_map can be removed I think but for now lets leave it
      if(i < ssa.size() - 1) {
        phi_x_inv[ssa[i+1].first] = i;
      }
    }

    // esa needs to be sorted for everything past this.
    std::sort(ssa.begin(), ssa.end()); // this wont be needed once sa_graph is sorted.
    std::sort(esa_sorted.begin(), esa_sorted.end());

    ulint i = 0; // ssa iterator
    ulint j = 0; // esa iterator
    ulint node = sa_map[esa_sorted.back()]; // init node as biggest value in pred ds because its circular.

    // computing the predecessor values using a merge-sort like combine phase
    while((i < ssa.size()) && (j < esa.size())) {
      if(ssa[i].first < esa_sorted[j]) {
        assert(node <= ssa[i].first);
        ulint successor_rank = pred.predecessor_rank_circular(phi_x_inv[ssa[i].first]) + 1;
        ulint successor = pred.select(successor_rank);
        // ulint predecessor = pred.select(successor_rank - 1);
        sa_graph[phi_x_inv[ssa[i].first]] = node; // this has to be improved, the hash is unnecessary
        bounds[phi_x_inv[ssa[i].first]].first = ssa[i].first - esa_sorted[j - 1]; // start_sample - pred(start_sample)
        bounds[phi_x_inv[ssa[i].first]].second = successor - esa[phi_x_inv[ssa[i].first]]; // successor of sample - pred of sample
        i += 1;
      }
      else {
        node = sa_map[esa_sorted[j]];
        j += 1;
      }
    }

    while(i < ssa.size()) {
      ulint successor_rank = pred.predecessor_rank_circular(esa[phi_x_inv[ssa[i].first]]) + 1;
      ulint successor = pred.select(successor_rank);
      ulint predecessor = pred.select(successor_rank - 1);
      sa_graph[phi_x_inv[ssa[i].first]] = node;
      bounds[phi_x_inv[ssa[i].first]].first = ssa[i].first - esa_sorted[j - 1];
      bounds[phi_x_inv[ssa[i].first]].second = successor - predecessor;
      i += 1;
    }
  }

  //! Build the csa/rads to be used with Phi.
  /*!
    \param unsorted_ssa Non-sorted start samples.
    \param ssa Sorted start samples.
    \param esa Non-sorted end samples.
    \param pred Predecessor bit vector. Used to calculate successors and predecessors.
  */
  void build_rads_phi(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    cout << "Building rads for phi ..." << endl;
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa; // we sort a temp array to not mess with esa
    sa_map.reserve(esa.size());
    sa_graph.resize(esa.size());
    bounds.resize(esa.size());

    cout << "Building SA map ..." << endl;
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[unsorted_ssa[i].first] = i; // sa_map just assigns run values to the ssa samples
      if(i > 0) {
        phi_x_inv[esa[i - 1]] = i; // array to tell you which samples phi takes you to
      }
    }

    ulint i = 0;
    ulint j = 0;
    ulint curr_pred = sa_map[ssa.back().first]; // assign pred to the last start sample
    std::sort(esa_sorted.begin(), esa_sorted.end()); // sort our temporary esa array

    // calculating predecessors and setting them as edges
    // sa_graph[...] = curr_pred; this means that the sample at the i-th run points to curr_pred
    cout << "Building graph ..." << endl;
    while((i < ssa.size()) && (j < esa.size())) {
      if(esa_sorted[j] < ssa[i].first) { // this will always be false on the first iteration. this means so long as the end sample is smaller than the start sample
        ulint sa_index = phi_x_inv[esa_sorted[j]];
        ulint successor = 0;
        ulint successor_rank = pred.predecessor_rank_circular(unsorted_ssa[sa_index].first) + 2;
        // ulint successor_rank = pred.predecessor_rank_circular(ssa[phi_x_inv[esa_sorted[i]]].first) + 1;
        // ulint predecessor = pred.select(successor_rank - 1);
        if(successor_rank >= pred.number_of_1()) { // this case happens when we are getting the successor of the largest sample
          successor = unsorted_ssa[sa_index].first + 1; // if true, then successor is just one bigger because the cost will be 1
        }
        else {
          successor = pred.select(successor_rank);
        }

        sa_graph[sa_index] = curr_pred; // sample -> curr_pred
        bounds[phi_x_inv[esa_sorted[j]]].first = esa_sorted[j] - ssa[i - 1].first; // lower_bound -> cost, sorted in SA order
        bounds[phi_x_inv[esa_sorted[j]]].second = successor - unsorted_ssa[sa_index].first; // upper_bound -> limit
        j += 1;
      }
      else {
        curr_pred = sa_map[ssa[i].first];
        i += 1;
      }
    }

    // if there are any leftover samples that need to point over to the last curr_pred that gets set, this is necessary
    while(j < esa.size()) {
      ulint sa_index = phi_x_inv[esa_sorted[j]];
      ulint successor = 0;
      ulint successor_rank = pred.predecessor_rank_circular(unsorted_ssa[sa_index].first) + 2;
      // ulint successor_rank = pred.predecessor_rank_circular(ssa[phi_x_inv[esa_sorted[i]]].first) + 1;
      // ulint predecessor = pred.select(successor_rank - 1);
      if(successor_rank >= pred.number_of_1()) {
        successor = unsorted_ssa[sa_index].first + 1;
      }
      else {
        successor = pred.select(successor_rank);
      }

      sa_graph[sa_index] = curr_pred;
      bounds[phi_x_inv[esa_sorted[j]]].first = esa_sorted[j] - ssa[i - 1].first;
      bounds[phi_x_inv[esa_sorted[j]]].second = successor - unsorted_ssa[sa_index].first;
      j += 1;
    }
  }

  //! Finds the cycles contained in the CSA graph.
  /*!
    \param ssa Start samples.
  */
  void find_cycles() {
    std::vector<uint> indegrees(sa_graph.size(), 0); // indegrees of the nodes
    std::vector<bool> visited(sa_graph.size(), false); // visited nodes so far
    auto temp_trees_bv = vector<bool>(sa_graph.size(), false);
    size_t tree_counter = 0;

    // counting indegrees of the nodes
    for(size_t i = 1; i < sa_graph.size(); i++) {
      if(sa_graph[i] >= 0)
        indegrees[sa_graph[i]] += 1;
    }

    // for all nodes with indegree 0, we check if they contain a cycle.
    for(size_t i = 1; i < sa_graph.size(); i++) {
      if(indegrees[i] != 0) { continue; }  //@ consider only starting nodes having no in-degree
      std::vector<ulint> current_path;
      ulint v = sa_graph[i];
      visited[i] = true;
      current_path.push_back(i);

      while(visited[v] == false && v != 0) {
        current_path.push_back(v);
        visited[v] = true;
        v = sa_graph[v];
      }

      // scan the current path and see if v is in it
      // bool is_cycle = false;
      // for(size_t path_it = 0; !is_cycle && path_it < current_path.size(); path_it++) {
      //   if(v == current_path[path_it]) {
      //     is_cycle = true;
      //   }
      // }

      // min threshold to remove really small cycles that arent worthwhile.

      // we can implement a min. path length threshold to include paths that
      // may not be cycles but can still be used to traverse samples.
      // if(is_cycle) { // if the path is a cycle we construct a tree
      if(current_path.size() >= 16) { // if the path has a length of at least 16 nodes, put it into a tree
        rads_tree branch = rads_tree<sparse_bv_type, rle_string_t>(current_path, bounds, tree_counter, tree_pointers);
        tree_counter += 1;
        trees.push_back(branch);
        for(size_t path_it = 0; path_it < current_path.size()-1; path_it++) { // -1 because the the last leaf node in our cycle is useless for queries
          temp_trees_bv[current_path[path_it]] = true; // set the nodes that are in the cycle to true in our bitvector
        }
      }
    }

    std::sort(tree_pointers.begin(), tree_pointers.end()); // sort based on the runs
    trees_bv = sparse_bv_type(temp_trees_bv);
    bounds.clear(); // clear because we dont need them anymore.
  }

  // sa_i is what we want to find
  // we can change std::vector<int> esa to int_vector<> esa (r_index copies samples_last_vec to samples_last)

  //! Query a particular SA[i] value.
  /*!
    \param sa_i The SA[i] to look for.
    \param bwt rle_string_t BWT needed.
    \param pred_to_run Bitvector used to find the run a particular predecessor is stored.
    \param sa Suffix array to query.
  */
  ulint query(ulint sa_i, rle_string_t &bwt, int_vector<> &pred_to_run, sparse_bv_type &pred, int_vector<> &sa, std::vector<ulint> &ssa) { // std::vector<ulint> &sa
    // if the sample we want is in a tree, use it, otherwise use phi
    // return SA value at position i (sa_i)
    const ulint initial_run = bwt.run_of_position(sa_i);
    const ulint run_l = bwt.run_range(initial_run).second; //@ the succeeding position we sampled
    ulint sa_j = sa[initial_run]; //@ = SA[initial_run], which we will update to SA[sa_i]
    ulint d = run_l-sa_i; //@ the distance of the closest sampled value up to the position `sa_i`

  //! Helper function that finds whether a sample is in a tree, therefore using it. If it is not in a tree it iterates Phi.
  /*!
    \param sa_j Sample at the end of the run that sa_i can be found in.
    \param d Distance between sa_j and sa_i.
    \param run Run in the r-index the sample is contained in.
    \param bwt BWT stored in the r-index.
    \param pred_to_run Bitvector providing us with the run of predecessors.
    \param pred Predecessor data structure from the r-index.
    \param sa Suffix array bitvector, in this case samples_last.
    \param ssa Unsorted start samples retrieved from the .ssa file post-construction.
  */
  // void helper_query(ulint &sa_j, ulint d, ulint run, rle_string_t &bwt, int_vector<> &pred_to_run, sparse_bv_type &pred, int_vector<> &sa, std::vector<ulint> &ssa) {
    ulint sa_jr;
    ulint sa_prime;
    ulint cost;
    ulint result;

    while(d > 0) {
      sa_jr = pred.predecessor_rank_circular(sa_j); // sample predecessor rank (run?)
      sa_prime = pred.select(sa_jr); // actual predecessor of sa_j
      cost = sa_j - sa_prime; // distance between sample and predecessor

      // check if pred is in a cycle
      // if on the last iteration (d == 1) just iterate phi
      if(trees_bv[pred_to_run[sa_jr]] && d != 1) {
        const auto current_run = pred_to_run[sa_jr];
        std::tuple<ulint, ulint, uint> tree_info = tree_pointers[trees_bv.rank(current_run)]; // get the tree that the sample belongs to.
        std::tuple<ulint, ulint, ulint> sa_prime_d_cost = trees[std::get<1>(tree_info)].query(std::get<2>(tree_info), cost, d); // tuple containing new sample run, distance travelled, and cost accumulated
        result = ssa[std::get<0>(sa_prime_d_cost)] + std::get<2>(sa_prime_d_cost); // sa_j is being set as the new sample
        sa_j = result;
        // run = pred_to_run[pred.predecessor_rank_circular(sa_j)]; //TODO: useless?
        d = d - std::get<1>(sa_prime_d_cost); // this is the distance left over.
      }
      else {
        // continue the iteration using phi
        ulint delta = sa_prime < sa_j ? sa_j - sa_prime : sa_j + 1; // distance between sample and its predecessor.
        ulint prev_sample = sa[pred_to_run[sa_jr] - 1]; // use samples_last (sa) to find the previous sample.
        sa_j = (prev_sample + delta) % bwt.size(); // get the next sample using delta and prev_sample.
        // run = pred_to_run[pred.predecessor_rank_circular(sa_j)]; // run of new sample. // TODO: useless?
        d -= 1;
      }
    }
  return sa_j;
  }

  void print_debug_info() {
    cout << "# of rads_trees & tree_pointers: " << trees.size() << endl;
    cout << "largest tree: " << get_largest_tree() << " samples" << endl;
  }

  void print_tree_runs(std::vector<std::pair<ulint, ulint>> &ssa, size_t start, size_t end) {
    for(size_t i = start; i < end; i++) {
      cout << "i: " << i << ", s_sample: " << ssa[std::get<0>(tree_pointers[i])].first << endl;
      cout << "run #: " << std::get<0>(tree_pointers[i]) << ", tree #: " << std::get<1>(tree_pointers[i]) << ", leaf node: " << std::get<2>(tree_pointers[i]) << endl;
    }
  }

  int get_largest_tree() {
    int max_tree_size = 0;
    for(int i = 0; i < trees.size(); i++) {
      if(trees[i].leaf_samples.size() > max_tree_size)
        max_tree_size = trees[i].leaf_samples.size();
    }

    return max_tree_size;
  }

  ulint serialize(std::ostream& out) {
    ulint w_bytes = 0;

    w_bytes += sdsl::serialize(trees.size(), out);
    for(size_t i = 0; i < trees.size(); i++) {
      w_bytes += trees[i].serialize(out);
    }

    w_bytes += sdsl::serialize(tree_pointers.size(), out);
    out.write((char*)tree_pointers.data(), tree_pointers.size()*sizeof(tree_pointers[0]));
    w_bytes += sizeof(tree_pointers[0])*tree_pointers.size();

    w_bytes += sdsl::serialize(bounds.size(), out);
    out.write((char*)bounds.data(), bounds.size()*sizeof(bounds[0]));
    w_bytes += sizeof(bounds[0])*bounds.size();

    w_bytes += sdsl::serialize(sa_graph.size(), out);
    out.write((char*)sa_graph.data(), sa_graph.size()*sizeof(sa_graph[0]));
    w_bytes += sizeof(sa_graph[0])*sa_graph.size();

    w_bytes += trees_bv.serialize(out);

    return w_bytes;
  }

  void load(std::istream& in) {
    size_t temp_size;

    in.read((char*)&temp_size, sizeof(temp_size));
    trees.resize(temp_size);
    for(size_t i = 0; i < temp_size; i++) {
      trees[i].load(in);
    }

    in.read((char*)&temp_size, sizeof(temp_size));
    tree_pointers.resize(temp_size);
    in.read((char*)tree_pointers.data(), temp_size*sizeof(tree_pointers[0]));

    in.read((char*)&temp_size, sizeof(temp_size));
    bounds.resize(temp_size);
    in.read((char*)bounds.data(), temp_size*sizeof(bounds[0]));

    in.read((char*)&temp_size, sizeof(temp_size));
    sa_graph.resize(temp_size);
    in.read((char*)sa_graph.data(), temp_size*sizeof(sa_graph[0]));

    trees_bv.load(in);
  }

protected:
  std::unordered_map<ulint, ulint> sa_map; // this is just pred_to_run but pred_to_run wasn't working?
  std::unordered_map<ulint, ulint> phi_x_inv; // map of phi or phi_inverse values

  // explain what tree pointers actually means better
  std::vector<rads_tree<>> trees; // list of trees.
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers; // pointers to the corresponding (run, tree, leaf node).
  std::vector<std::pair<ulint,ulint>> bounds; // lower and upper bounds of each node in the sa graph. // can be deleted at some point
  std::vector<ulint> sa_graph; // adj. list representing the sa graph.
  sparse_bv_type trees_bv; // bitvector that tells us which samples are in trees.
};
}

#endif
