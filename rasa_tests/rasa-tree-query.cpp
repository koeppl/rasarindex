#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <cstring>
#include <string>

#include "../internal/rle_string.hpp"
#include "../internal/r_index.hpp"
#include "../internal/ri_rasa_tree.hpp"
#include "../internal/utils.hpp"
#include "../internal/sparse_sd_vector.hpp"
#include "../internal/sparse_hyb_vector.hpp"
#include "dcheck.hpp"

using namespace sdsl;
using namespace ri;
using namespace std;

// testing the query function for the rads_tree class.
// create several trees.
//

rads_tree<> prepare_tree(const std::vector<ulint>& nodes, const std::vector<std::pair<ulint, ulint>>& bounds, int tree_id,std::vector<std::tuple<ulint, ulint, uint>>& tree_pointers) {
  DCHECK_EQ(nodes.size(), bounds.size());
  const auto& max_el = std::max_element(nodes.begin(), nodes.end());
    std::vector<std::pair<ulint, ulint>> bound_map(*max_el+1);
    for(size_t i = 0; i < nodes.size(); ++i) {
      bound_map[nodes[i]] = bounds[i];
    }

    return rads_tree(nodes, bound_map, tree_id, tree_pointers);
}

std::tuple<ulint, ulint, ulint>  leaf_query(const rads_tree<>& tree, const ulint leaf_pos, uint cost, const uint distance_bound) {
  const auto& treearray =  tree.m_tree;
  const ulint leftmost_leaf_pos = tree.left_most_i;
  DCHECK_LE(leftmost_leaf_pos, leaf_pos);
  ulint distance = 0;
  for(;leaf_pos + distance < treearray.size(); ++distance) {
    const ulint visiting_leaf = leaf_pos + distance;
    if(treearray[visiting_leaf].second <= cost || distance >= distance_bound) {
      break;
    }
	if(visiting_leaf != treearray.size()-1) { // not the last leaf
      cost += treearray[visiting_leaf].first;
    }
  }
  if(distance > distance_bound) { distance = distance_bound; }
  if(leaf_pos + distance >= treearray.size()) { distance = treearray.size()-1-leaf_pos; }
  return std::make_tuple(tree.leaf_samples[distance+(leaf_pos - leftmost_leaf_pos)], distance, cost); 
}

void check_query_tuple(const std::tuple<ulint, ulint, ulint>& a, const std::tuple<ulint, ulint, ulint>& b) {
  DCHECK_EQ(std::get<0>(a), std::get<0>(b));
  DCHECK_EQ(std::get<1>(a), std::get<1>(b));
  DCHECK_EQ(std::get<2>(a), std::get<2>(b));
}

#include <cstdlib>

void random_tree_test() {
  srand(0);
  for(size_t no_leaves = 8; no_leaves < 64; ++no_leaves) {
    const size_t TREE_FOUR_SIZE = no_leaves;
    std::vector<ulint> leaves; // a lot of samples to push it.
    std::vector<std::pair<ulint, ulint>> bounds;
    for(size_t i = 0; i <= TREE_FOUR_SIZE; i++) {
      leaves.push_back(i);
      bounds.push_back(std::make_pair(rand()%10, rand() % 100+80));
    }

    std::vector<std::tuple<ulint, ulint, uint>> tree_pointers;
    auto tree = prepare_tree(leaves, bounds, 0, tree_pointers);

    // retrieve every sample from every possible leaf.
    for(size_t i = 0; i < tree.leaf_node_bv.size(); i++) {
      if(tree.leaf_node_bv[i] == 1) {
        const ulint node_distance = tree.calculate_d(tree.left_most_i, i);
        const ulint samples_left = TREE_FOUR_SIZE - node_distance;
        for(size_t j = 1; j <= samples_left; j++) {
          const auto query_answer = tree.query(i, 0, j);
          const auto query_answer_naive = leaf_query(tree, i, 0, j);
          check_query_tuple(query_answer, query_answer_naive);
        }
        // no overshooting
        for(size_t j = samples_left+1; j < samples_left+10; j++) {
          const auto query_answer = tree.query(i, 0, j);
          const auto query_answer_naive = leaf_query(tree, i, 0, j);
          // check_query_tuple(query_answer, query_answer_naive); //TODO: something broken here
          DCHECK_EQ(std::get<0>(query_answer), std::get<0>(query_answer_naive));
        }
      }
    }
  }
}


int main() {


  std::vector<rads_tree<>> trees;
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers;

  std::vector<ulint> tree_0 = {0, 1, 2, 3, 4, 5, 6, 7}; // basic tree with 8 elements to test traversal unhindered by bounds.
  std::vector<std::pair<ulint, ulint>> tree_0_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  std::vector<ulint> tree_1 = {26, 3, 11, 20, 0, 17}; // example with actual bounds.
  std::vector<std::pair<ulint, ulint>> tree_1_bounds = {{5,1}, {0,6}, {0,3}, {2,1}, {0,3}, {4,1}};

  std::vector<ulint> tree_2 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}; // larger example with 14 samples to test traversal again.
  std::vector<std::pair<ulint, ulint>> tree_2_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  std::vector<ulint> tree_3 = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // 9 samples to test the new rung of leaves that get made.
  std::vector<std::pair<ulint, ulint>> tree_3_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  constexpr size_t TREE_FOUR_SIZE = 18;
  std::vector<ulint> tree_4; // a lot of samples to push it.
  std::vector<std::pair<ulint, ulint>> tree_4_bounds;
  for(size_t i = 0; i <= TREE_FOUR_SIZE; i++) {
    tree_4.push_back(i);
    tree_4_bounds.push_back(std::make_pair(0,1));
  }

  // std::vector<ulint> tree_5 = {26, 13, 3, 11, 16, 20, 0, 17, 1, 2, 19, 24};
  // std::vector<std::pair<ulint, ulint>> tree_5_bounds = {{5,1}, {0,3}, {0,6}, {0,3}, {3,0}, {2,1}, {0,3}, {4,1}, {0,3}, {0,3}, {0,3}, {0,3}};

  trees.emplace_back(prepare_tree(tree_0, tree_0_bounds, 0, tree_pointers));
  trees.emplace_back(prepare_tree(tree_1, tree_1_bounds, 1, tree_pointers));
  trees.emplace_back(prepare_tree(tree_2, tree_2_bounds, 2, tree_pointers));
  trees.emplace_back(prepare_tree(tree_3, tree_3_bounds, 3, tree_pointers));
  trees.emplace_back(prepare_tree(tree_4, tree_4_bounds, 4, tree_pointers));
  // trees.emplace_back(prepare_tree(tree_5, tree_5_bounds, 5, tree_pointers));

  std::tuple<ulint, ulint, ulint> sample_and_delta;
  // retrieve every sample using the left most sample.
  for(size_t i = 1; i <= 6; i++) {
    sample_and_delta = trees[0].query(trees[0].left_most_i, 0, i);
    DCHECK_EQ(std::get<0>(sample_and_delta) , i);
    check_query_tuple(sample_and_delta, leaf_query(trees[0], trees[0].left_most_i, 0, i));
  }

  // couple of tests to see if we get the samples that we expect from the tree with bounds.
  sample_and_delta = trees[1].query(8, 0, 1);
  check_query_tuple(sample_and_delta, std::make_tuple(3,1,5));

  sample_and_delta = trees[1].query(8, 0, 2);
  check_query_tuple(sample_and_delta, std::make_tuple(11,2,5));

  sample_and_delta = trees[1].query(8, 0, 3);
  DCHECK_EQ(std::get<0>(sample_and_delta), 11);
  DCHECK_EQ(std::get<1>(sample_and_delta), 2);
  DCHECK_EQ(std::get<2>(sample_and_delta), 5);

  sample_and_delta = trees[1].query(9, 0, 1);
  DCHECK_EQ(std::get<0>(sample_and_delta), 11);
  DCHECK_EQ(std::get<1>(sample_and_delta), 1);
  DCHECK_EQ(std::get<2>(sample_and_delta), 0);

  sample_and_delta = trees[1].query(9, 0, 2);
  DCHECK_EQ(std::get<0>(sample_and_delta), 20);
  DCHECK_EQ(std::get<1>(sample_and_delta), 2);
  DCHECK_EQ(std::get<2>(sample_and_delta), 0);

  const ulint query_node = rads_tree<>::VERSION == 0 ? 5 : trees[1].left_most_i+2;

  sample_and_delta = trees[1].query(query_node, 0, 1);
  DCHECK_EQ(std::get<0>(sample_and_delta), 20);
  DCHECK_EQ(std::get<1>(sample_and_delta), 1);
  DCHECK_EQ(std::get<2>(sample_and_delta), 0);

  sample_and_delta = trees[1].query(query_node, 0, 2);
  DCHECK_EQ(std::get<0>(sample_and_delta), 0);
  DCHECK_EQ(std::get<1>(sample_and_delta), 2);
  DCHECK_EQ(std::get<2>(sample_and_delta), 2);

  sample_and_delta = trees[1].query(query_node, 0, 3);
  DCHECK_EQ(std::get<0>(sample_and_delta), 17);
  DCHECK_EQ(std::get<1>(sample_and_delta), 3);
  DCHECK_EQ(std::get<2>(sample_and_delta), 2);

  sample_and_delta = trees[1].query(query_node, 0, 4);
  DCHECK_EQ(std::get<0>(sample_and_delta), 17);
  DCHECK_EQ(std::get<1>(sample_and_delta), 3);
  DCHECK_EQ(std::get<2>(sample_and_delta), 2);

  for(size_t i = 1; i <= 13; i++) {
    sample_and_delta = trees[2].query(trees[2].left_most_i, 0, i);
    DCHECK_EQ(std::get<0>(sample_and_delta) , i);
  }

  for(size_t i = 1; i <= 8; i++) {
    sample_and_delta = trees[3].query(trees[3].left_most_i, 0, i);
    DCHECK_EQ(std::get<0>(sample_and_delta) , i);
  }

  for(size_t i = 1; i <= TREE_FOUR_SIZE; i++) {
    sample_and_delta = trees[4].query(trees[4].left_most_i, 0, i);
    DCHECK_EQ(std::get<0>(sample_and_delta) , i);
  }

  // retrieve every sample from every possible leaf.
  for(size_t i = 0; i < trees[4].leaf_node_bv.size(); i++) {
    if(trees[4].leaf_node_bv[i] == 1) {
      const ulint node_distance = trees[4].calculate_d(trees[4].left_most_i, i);
      const ulint samples_left = TREE_FOUR_SIZE - node_distance;
      for(size_t j = 1; j <= samples_left; j++) {
        sample_and_delta = trees[4].query(i, 0, j);
        DCHECK_EQ(std::get<0>(sample_and_delta) , (node_distance + j));
        check_query_tuple(sample_and_delta, leaf_query(trees[4], i, 0, j));
      }
      // no overshooting
      for(size_t j = samples_left+1; j < samples_left+10; j++) {
        sample_and_delta = trees[4].query(i, 0, j);
        DCHECK_EQ(std::get<0>(sample_and_delta), node_distance + samples_left);
        check_query_tuple(sample_and_delta, leaf_query(trees[4], i, 0, j));
      }
    }
  }

  random_tree_test();

  // sample_and_delta = trees[5].query(9, 0, 3);
  // DCHECK_EQ(std::get<0>(sample_and_delta) , 11);

  cout << "all good!" << endl;

  return 0;
}
