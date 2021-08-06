#ifndef R_INDEX_R_A_
#define R_INDEX_R_A_

#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <definitions.hpp>
#include <r_index.hpp>
#include <ri_rasa_tree.hpp>

using namespace sdsl;

namespace ri {
class rads {
public:
  /*
   * get comments in on everything
   * with comments you get the pseudocode written, the paper written!
   * easier to revisit
   */

  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    init_by_value(ssa, esa);
    cout << "Done. Now listing paths." << endl;
    list_paths(ssa, esa);
  }

  void init_by_value(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    assert(ssa.size() == esa.size());
    esa_map.reserve(esa.size());
    phi_inv_sa.resize(esa.size());
    bounds.resize(esa.size());
    trees_sample_bv.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      esa_map[esa[i]] = i;
      if(i < ssa.size() - 1) {
        pis_inv[ssa[i+1].first] = i;
      }
    }

    std::sort(ssa.begin(), ssa.end()); // this wont be needed once phi_inv_sa is sorted.
    std::sort(esa.begin(), esa.end());

    ulint i = 0; // ssa iterator
    ulint j = 0; // esa iterator
    ulint node = esa_map[esa.back()]; // init node as biggest value in pred ds because its circular.

    // computing the predecessor values using a merge-sort like combine phase
    while((i < ssa.size()) && (j < esa.size())) {
      if(ssa[i].first < esa[j]) {
        assert(node <= ssa[i].first);
        phi_inv_sa[pis_inv[ssa[i].first]] = node; // this has to be improved, the hash is unnecessary
        bounds[pis_inv[ssa[i].first]].first = ssa[i].first - esa[j-1];
        bounds[pis_inv[ssa[i].first]].second = esa[j] - ssa[i].first;
        i += 1;
      }
      else {
        node = esa_map[esa[j]];
        j += 1;
      }
    }

    while(i < ssa.size()) {
      phi_inv_sa[pis_inv[ssa[i].first]] = node;
      bounds[pis_inv[ssa[i].first]].first = ssa[i].first - esa[j-1];
      bounds[pis_inv[ssa[i].first]].second = esa[j] - ssa[i].first;
      i += 1;
    }
  }

  // list_paths() finds all cycles in our graph.
  void list_paths(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    std::vector<uint> indegrees(phi_inv_sa.size(), 0); // indegrees of the nodes
    std::vector<bool> visited(phi_inv_sa.size(), false); // visited nodes so far

    // counting indegrees of the nodes
    for(size_t i = 0; i < phi_inv_sa.size(); i++) {
      if(phi_inv_sa[i] >= 0)
        indegrees[phi_inv_sa[i]] += 1;
    }

    // for all nodes with indegree 0, we check if they are a cycle.
    for(size_t i = 0; i < phi_inv_sa.size(); i++) {
      if(indegrees[i] == 0) {
        std::vector<uint> current_path;
        int u = i;
        int v = phi_inv_sa[u];
        visited[u] = true;
        current_path.push_back(u);

        while(visited[v] == false) {
          current_path.push_back(v);
          visited[v] = true;
          v = phi_inv_sa[v];
        }

        // scan the current path and see if v is in it
        bool is_cycle = false;
        for (size_t i = 0; !is_cycle && i < current_path.size(); i++) {
          if(v == current_path[i])
            is_cycle = true;
        }

        // we can implement a min. path length threshold to include paths, that
        // may not be cycles but can still be used to traverse samples.
        if(is_cycle) { // if the path is a cycle we construct a tree
          rads_tree branch = rads_tree(current_path, bounds);
          trees.push_back(branch);
          for(size_t i = 0; i < current_path.size(); i++) {
            trees_sample_bv[current_path[i]] = true; // set the nodes that are in the cycle to true in our bitvector
          }
        }
      }
    }
  }

  inline int get_size() {
    return phi_inv_sa.size();
  }

  inline int get_num_paths() {
    return trees.size();
  }

protected:
  std::unordered_map<ulint, ulint> esa_map;
  std::unordered_map<ulint, ulint> pis_inv;

  std::vector<std::pair<ulint,ulint>> bounds; // lower and upper bounds of each node in the sa graph.
  std::vector<ulint> phi_inv_sa; // adj. list representing the sa graph.
  std::vector<rads_tree> trees; // list of cycle trees.
  std::vector<bool> trees_sample_bv; // bitvector that tells us which samples are in trees.
};
}

#endif
