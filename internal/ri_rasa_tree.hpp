#ifndef R_INDEX_R_A_TREE_
#define R_INDEX_R_A_TREE_

#include <definitions.hpp>
#include <r_index.hpp>
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

using namespace sdsl;
namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>

class rads_tree {
public:
  std::vector<std::pair<long long int, long long int>> tree; // nodes are pairs that represent: (edge cost, edge threshold).
  std::vector<ulint> leaf_samples; // actual sa sample values represented by run index
  sparse_bv_type leaf_node_bv; // bv telling us which node is a leaf node
  uint left_most_i;

  rads_tree(){};
  rads_tree(std::vector<ulint> &cycle, std::vector<std::pair<ulint, ulint>> &bounds, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers){
    size_t n = cycle.size(); // path size
    leaf_samples = cycle; // store the samples in the cycle
    tree = std::vector<std::pair<long long int, long long int>>((size_t)1<<(size_t)(ceil(log2(n)) + 1), std::make_pair(-1,0)); // tree size initialization, default the costs to -1
    auto temp_leaf_bv = vector<bool>(tree.size(), false);

    constructor_helper(bounds, 1, 0, n-1, tree_num, tree_pointers, temp_leaf_bv);
    leaf_node_bv = sparse_bv_type(temp_leaf_bv);
    left_most_i = (tree.size()>>1);
  }

  // so many arguments...
  // recursively constructs a balanced binary tree from the leaves up to the root.
  void constructor_helper(std::vector<std::pair<ulint, ulint>> &bounds, size_t node, size_t begin, size_t end, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers, std::vector<bool> &leaf_bv) {
    size_t mid = (begin + end)/2;
    if(begin >= end) { // condition which creates the leaf nodes.
      tree[node].first = (long long int)bounds[begin].first; // edge cost
      tree[node].second = (long long int)bounds[begin].second; // edge threshold
      tree_pointers.push_back(std::make_tuple(leaf_samples[begin], tree_num, node));
      //cout << leaf_samples[begin] << " " << node << endl;
      leaf_bv[node] = true;
      return;
    }

    constructor_helper(bounds, 2*node, begin, mid, tree_num, tree_pointers, leaf_bv); // construct left child
    constructor_helper(bounds, 2*node+1, mid+1, end, tree_num, tree_pointers, leaf_bv); // construct right child
    tree[node].first = tree[2*node].first + tree[2*node+1].first; // edge sum of left and right child
    tree[node].second = std::min(tree[2*node].second, tree[2*node+1].second - tree[2*node].first); // edge threshold
  }

  // args: leaf_pos, cost, d | returns: returns new sample and distance left
  //                           (sa', d)
  std::pair<ulint, ulint> query(ulint leaf_pos, uint cost, uint d) {
    // begin climbing the tree and only collect when we go to a sibling.
    // when you go to a parent do not collect the money, that is if the parent has the same left most node.
    // if youre moving to a subtree where the leftmost node is not the same then we have to collect the money.
    cout << "Querying tree ..." << endl;
    cout << "leaf_pos: " << leaf_pos << ", cost: " << cost << ", d: " << d << endl;
    return (climb(leaf_pos, cost, d));
  }

  // climbing conditions
  // 1. node is non-negative
  // 2. cost does not exceed nodes threshold
  std::pair<ulint, ulint> climb(ulint &node_pos, uint &cost, uint d) {
    // this provides us with the index of the sample that is being used to query
    // left_most_i is 0 and we get the distance from our start_pos
    // at the end we do leaf_index + the distance that was travelled to get the new sample
    // eg. start_pos = 21 -> leaf_index = 0; distance travelled = 2; leaf_index + distance travelled = 2; new_sample = leaf_samples[2] (11)
    cout << "Climbing ..." << endl;
    cout << "!pre calc_d!" << endl;
    ulint start_pos = node_pos;
    ulint start_leaf_index = calculate_d(left_most_i, start_pos); // this isnt necessary because well we know the index we start at right but nice way to test calc_d for now
    cout << "!post calc_d!" << endl;

    // check if when we shift up we get to the root of the left subtree
    // meaning we can skip forward to the beginning of the right subtree and descend
    if(node_pos>>(__builtin_ctzl(~node_pos)) == 2) {
      cout << "\nSkip to the right subtree ..." << endl;
      node_pos = 3;

      //cost += tree[node_pos-1].first; // add cost of node that denies us | ? adding this cost is unnecessary ?

      descend(node_pos, cost);
      ulint distance = calculate_d(start_pos, node_pos);

      cout << "leaf sample: " << leaf_samples[start_leaf_index + distance] + cost << endl;
      cout << "distance left: " << d - distance << endl;

      return (std::make_pair(leaf_samples[start_leaf_index + distance] + cost, d - distance)); // return new sample and new distance
    }

    cout << "!pre climb!" << endl;
    while((node_pos != 1) && (tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) { // climb up while the upper_bounds let us
      node_pos = node_pos >> 1; 
    }
    cout << "!post climb!" << endl;

    cout << "we're at the root or we've been denied." << endl;
    cost += tree[node_pos<<1].first; // add cost of node that denies us
    node_pos = (node_pos<<1) + 1;

    descend(node_pos, cost); // descend using new node_pos, cost, and d
    ulint distance = calculate_d(start_pos, node_pos);

    cout << "leaf sample: " << leaf_samples[start_leaf_index + distance] + cost << endl;
    cout << "distance left: " << d - distance << endl;

    return (std::make_pair(leaf_samples[start_leaf_index + distance] + cost, d - distance)); // return new sample and new distance
  }

  // during the descent the node_pos gets shifted around so that calculate_d can do its job
  void descend(ulint &node_pos, uint &cost) {
    // 1. check node_pos that is passed in
    // 2. if ✓ then check your immediate sibling
    // 3. if ✕ then check left child until you get a green light
    // while((node_pos < tree.size()) && (tree[node_pos].first >= 0)) { // 1. check if we are out of bounds 2. check if there is even a node here
    cout << "Descending ..." << endl;
    while((node_pos < leaf_node_bv.size()) && !leaf_node_bv[node_pos]) { // while we are not at a leaf node
      if(node_pos & 1 == 0) { // this means we are at a left child node
        cout << "left node." << endl;
        if((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) { // if good -> go to right sibling
          cout << "go to right sibling." << endl;
          cost += tree[node_pos].first;
          node_pos += 1;
        }
        else {
          cout << "go to left child." << endl;
          node_pos = node_pos << 1; // if not good then go to left child and check on next iteration
        }
      }
      else { // right node
        cout << "right node." << endl;
        if((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) { // if good -> go to right child
          cout << "go to right child." << endl;
          node_pos = (node_pos << 1) + 1;
        }
        else {
          cout << "go to left child." << endl;
          node_pos = node_pos << 1; // if not good then go to left child and check that guy on the next iterations
        }
      }
    }

    // cout << "node pos before shift: " << node_pos << endl;
    // node_pos = node_pos >> 1;
    cout << "node pos after shifts: " << node_pos << endl;
  }

  // calculates distance between leaf nodes start_pos and end_pos
  ulint calculate_d(ulint start_pos, ulint end_pos) {
    cout << "\nCalculate distance between nodes ..." << endl;
    cout << "start pos: " << start_pos << ", end pos: " << end_pos << endl;
    ulint d = 0;

    // cout << "\nBit vector:" << endl;
    // for(int i = 0; i < leaf_node_bv.size(); i++) {
    //   cout << leaf_node_bv[i] << endl;
    // }

    if(start_pos >= left_most_i) { // starting from bottom layer
      cout << "start pos >= left_most_i." << endl;
      if(end_pos >= left_most_i) { // going to bottom layer bot->bot
        cout << "end_pos >= left_most_i" << endl;
        cout << "end pos rank" << endl;
        d += leaf_node_bv.rank(end_pos);
        cout << "start pos rank" << endl;
        d -= leaf_node_bv.rank(start_pos);
        cout << "start pos shift rank" << endl;
        d -= leaf_node_bv.rank((start_pos>>1));
        cout << "end pos shift rank" << endl;
        d += leaf_node_bv.rank((end_pos>>1));
      }
      else { // going to top layer bot->top
        cout << "endpos < left_most_i" << endl;
        d += leaf_node_bv.rank(end_pos) - leaf_node_bv.rank(start_pos) - leaf_node_bv.rank((start_pos>>1)) + leaf_node_bv.rank((end_pos<<1));
      }
    }
    else { // starting from top layer
      cout << "start pos < left_most_i" << endl;
      if(end_pos >= left_most_i) {// going to bottom layer top->bot
        cout << "end pos >= left_most_i" << endl;
        d += leaf_node_bv.rank(end_pos) - leaf_node_bv.rank(start_pos) - leaf_node_bv.rank((start_pos<<1)) + leaf_node_bv.rank((end_pos>>1));
      }
      else { // going to top layer top->top
        cout << "end pos < left_most_i" << endl;
        d += leaf_node_bv.rank(end_pos) - leaf_node_bv.rank(start_pos) - leaf_node_bv.rank((start_pos<<1)) + leaf_node_bv.rank((end_pos<<1));
      }
    }

    cout << "\nd: " << d << endl;
    return d;
  }

  void print_array() {
    cout << "Printing array ..." << endl;
    for (size_t i = 0; i < tree.size(); i++) {
      cout << "(" << tree[i].first << ", " << tree[i].second << ") ";
    }

    cout << endl;
    int leaf_counter = 0;
    for (size_t i = 0; i < leaf_node_bv.size(); i++) {
      if(leaf_node_bv[i]) {
        cout << " 1,";
        leaf_counter++;
      }
      else {
        cout << " -,";
      }
    }
    cout << endl;
  }
};
}

#endif
