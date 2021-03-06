#ifndef CSATREE_CLASSIC_HPP
#define CSATREE_CLASSIC_HPP

#include "definitions.hpp"
#include "r_index.hpp"
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

using namespace sdsl;

namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>

//! Tree data structure that stores the samples contained within a cycle. Performs queries.
class CSATree {
public:
  using cost_type = ulint;
  using limit_type = long long int;

  std::vector<std::pair<cost_type, limit_type>> tree; // nodes are pairs that represent: (edge cost, edge threshold).
  std::vector<ulint> leaf_samples; // actual sa sample values represented by run index
  sparse_bv_type leaf_node_bv; // bv telling us which node is a leaf node
  uint left_most_i; // index at which the left most leaf is stored || sdsl::serialize
  size_t height; // height of the tree starting at 0 || sdsl::serialize

  CSATree(){};
  //! CSATree constructor that recursively builds a balanced binary interval tree using the samples provided.
  /*!
    \param cycle Path of sample nodes on which to be built a tree.
    \param bounds Calculated bounds used for lower and upper bound calculation.
    \param tree_num Number of tree.
    \param tree_pointers Vector of tree pointers that is added to every time a leaf is set.
  */
  CSATree(const std::vector<ulint> &cycle, const std::vector<std::pair<ulint, ulint>> &bounds, const uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers) {
    size_t n = cycle.size(); // path size
    height = ceil(log2(n)); // height of the tree
    leaf_samples = cycle; // store the samples in the cycle
    tree = decltype(tree)(1ULL<<(height + 1), std::make_pair(-1, 0)); // tree size initialization, default the costs to -1
    auto temp_leaf_bv = vector<bool>(tree.size(), false);

    constructor_helper(bounds, 1, 0, n - 1, tree_num, tree_pointers, temp_leaf_bv);
    leaf_node_bv = sparse_bv_type(temp_leaf_bv);
    left_most_i = (tree.size() >> 1);
  }

  CSATree(const CSATree &other_tree) {
    this->tree = other_tree.tree;
    this->leaf_samples = other_tree.leaf_samples;
    this->leaf_node_bv = other_tree.leaf_node_bv;
    this->left_most_i = other_tree.left_most_i;
    this->height = other_tree.height;
  }

  CSATree(CSATree &&other_tree) noexcept
  : tree(move(other_tree.tree))
  , leaf_samples(move(other_tree.leaf_samples))
  , leaf_node_bv(move(other_tree.leaf_node_bv))
  , left_most_i(move(other_tree.left_most_i))
  , height(move(other_tree.height))
  {}

  CSATree& operator=(const CSATree &other_tree) {
    return *this = CSATree(other_tree);
  }

  CSATree& operator=(CSATree &&other_tree) {
    swap(tree, other_tree.tree);
    swap(leaf_samples, other_tree.leaf_samples);
    swap(leaf_node_bv, other_tree.leaf_node_bv);
    swap(left_most_i, other_tree.left_most_i);
    swap(height, other_tree.height);
    return *this;
  }

  //! Helper to the constructor that goes to left and right children and sets leaf nodes once they have been encountered, and calculates the parents cost & bounds.
  /*!
    \param bounds Sample bounds that were calculated during construction.
    \param node Tree node that we are currently on.
    \param begin Left value that we are covering in the current interval.
    \param end Right value that we are covering in the current interval.
    \param tree_num Tree # we are building currently.
    \param tree_pointers Pointers that are added to for eveyr leaf we set.
    \param leaf_bv Bitvector to know which nodes are leaf nodes.
  */
  void constructor_helper(const std::vector<std::pair<ulint, ulint>> &bounds, const size_t node, const size_t begin, const size_t end, const uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers, std::vector<bool> &leaf_bv) {
    const size_t mid = (begin + end)/2;
    if(begin >= end) { // condition which creates the leaf nodes.
      tree[node].first = bounds[leaf_samples[begin]].first; // edge cost
      tree[node].second = bounds[leaf_samples[begin]].second; // edge threshold
      if(begin != leaf_samples.size()-1) { // if it is not the last sample in the tree, add a tree pointer to be used when querying
        tree_pointers.push_back(std::make_tuple(leaf_samples[begin], tree_num, node));
      }
      leaf_bv[node] = true;
      return;
    }

    constructor_helper(bounds, 2*node, begin, mid, tree_num, tree_pointers, leaf_bv); // construct left child
    constructor_helper(bounds, 2*node+1, mid+1, end, tree_num, tree_pointers, leaf_bv); // construct right child
    tree[node].first = tree[2*node].first + tree[2*node+1].first; // edge sum of left and right child
    tree[node].second = std::min<limit_type>(tree[2*node].second, tree[2*node+1].second - tree[2*node].first); // edge threshold
  }

  //! Query function.
  /*!
    \param node_pos Leaf node to begin traversal.
    \param cost Cost that we will carry and add to as we climb.
    \param d Distance that we would ideally travel.
  */
  std::tuple<ulint, ulint, ulint> query(ulint node_pos, uint cost, uint d) {
    //! Climbing traversal function. Climb until we cannot and then begin descending.
    // begin climbing the tree and only collect when we go to a sibling.
    // when you go to a parent do not collect the cost, that is if the parent has the same left most node.
    // if youre moving to a subtree where the leftmost node is not the same then we have to collect the money.
   
    ulint current_height = height;
    ulint start_pos = node_pos;
    ulint prev_pos = node_pos;
    ulint max_d_travelled = 0;

    // we can enter at either the height or height - 1
    if(node_pos < left_most_i) {
      current_height -= 1;
    }

    // climb up while the upper_bounds let us / climb up at least once if we are still at the start node
    ulint current_leftmost_node = (node_pos << (height - current_height));
    while(((node_pos != 1) && (tree[node_pos].second > 0) && (cost < tree[node_pos].second)) || (start_pos == node_pos)) {
      const auto last_bit = (node_pos & 1);

      // if we are on the right branch of nodes do not climb anymore
      // we also check if the sample we want is in the current interval
      if((d <= max_d_travelled) || (static_cast<ulint>(__builtin_ctzl(~node_pos)) == (current_height + 1))) { // find lowest ancestor for which we are in the left subtree. All ancestors between do not give us more information.
        descend(start_pos, node_pos, cost, d, current_height);
        ulint distance = calculate_d(start_pos, node_pos);
        ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
        return (std::make_tuple(leaf_samples[leaf_sample_index], distance, cost)); // return new sample and new distance
      }

      // check if when we shift up we get to the root of the left subtree
      // meaning we're on the rightmost branch of the left subtree
      // this means we can skip forward to the right subtree
      prev_pos = node_pos;
      if(node_pos>>(__builtin_ctzl(~node_pos)) == 2) {
        node_pos = 3;
        current_height = 1;
        cost += tree[prev_pos].first;
        descend(start_pos, node_pos, cost, d, current_height);
        ulint distance = calculate_d(start_pos, node_pos);
        ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
        return (std::make_tuple(leaf_samples[leaf_sample_index], distance, cost)); // return new sample and new distance
      }

      // left node up
      if(last_bit == 0) {
        node_pos = node_pos >> 1;
      }
      else {
        // right node goes up and right
        if(node_pos != start_pos) { // check parent node if we are no longer at our start position
          if((tree[node_pos >> 1].second > 0) && (cost < tree[node_pos >> 1].second)) {
            node_pos = (node_pos >> 1) + 1;
          }
          else {
            break; // get out cause you're going to make an illegal jump
          }
        }
        else {
          node_pos = (node_pos >> 1) + 1;
        }
      }

      current_height -= 1;

      // moved to a subtree where the left most node is not the same, therefore we collect
      if((node_pos << (height - current_height)) != current_leftmost_node) {
        cost += tree[prev_pos].first;
        current_leftmost_node = (node_pos << (height - current_height));
      }

      const ulint max_node_pos = ((node_pos << (height - current_height)) + ((1 << (height - current_height)) - 1)); // furthest right node we can reach from current position
      max_d_travelled = calculate_d(start_pos, max_node_pos); // (node_pos gets shifted by the distance to the leaves) + (1 shifted that many times left - 1)
    }

    // if we climbed up from a right node, we might encounter a node that denies us,
    // and the nodes below that might also deny us, therefore we shouldn't collect.
    // this descends until we have found a node that lets us pass and we collect that.
    const auto last_bit = (prev_pos & 1);
    if(last_bit == 1) {
      while(!leaf_node_bv[node_pos]) {
        if(!(tree[node_pos].second > 0) || !(cost < tree[node_pos].second)) { // this means we keep getting denied
          node_pos = (node_pos << 1);
          current_height += 1;
        }
        else {
          node_pos = (node_pos << 1) + 1; // once we find a bound that were allowed to go through, go to the right child.
          current_height += 1;
          cost += tree[node_pos - 1].first;
          descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
          break;
        }
      }
    }
    else {
      cost += tree[node_pos << 1].first; // add cost of node that denies us
      node_pos = (node_pos << 1) + 1; // move node_pos to the right child of the current node // or should it be the left child of the node
      current_height += 1;
      descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
    }

    ulint distance = calculate_d(start_pos, node_pos); // calculate the distance from start to end to know how far we travelled
    ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
    return (std::make_tuple(leaf_samples[leaf_sample_index], distance, cost)); // return new sample and new distance
  }

  //! Descent function. This descends until we reach a leaf node making the proper movements based on our cost and the bounds we encounter.
  /*!
    \param start_pos Node we started the whole traversal from.
    \param node_pos Node position that we are currently at.
    \param cost Cost we have accumulated so far.
    \param current_height Current height in the tree that we're at.
  */
  void descend(ulint start_pos, ulint &node_pos, uint &cost, uint d, ulint current_height) {
    ulint prev_pos = node_pos;

    while((node_pos < leaf_node_bv.size()) && !leaf_node_bv[node_pos]) { // while we are in the bounds of the tree and not at a leaf node
      const auto last_bit = (node_pos & 1);
      prev_pos = node_pos;
      if(last_bit == 0) { // left node descent
        if((tree[node_pos].second > 0) && (cost < tree[node_pos].second)) { // if good, go to right sibling
          ulint min_node = ((node_pos + 1) << (height - current_height));
          const auto min_d_travelled = calculate_d(start_pos, min_node);
          if(d < min_d_travelled) {
            // this means we would go too far when going to right sibling
            // go to right child instead
            node_pos = node_pos << 1;
            node_pos += 1;
            current_height += 1;
            cost += tree[node_pos - 1].first;
          }
          else {
            // go to right sibling
            node_pos += 1;
            cost += tree[prev_pos].first;
          }
        }
        else {
          // if not good then go to left child and check on next iteration
          node_pos = node_pos << 1;
          current_height += 1;
        }
      }
      else { // right node
        if((tree[node_pos].second > 0) && (cost < tree[node_pos].second)) { // if good, go to right child
          ulint min_node = (node_pos << (height - current_height)); // this is how you can get the left mode node of the current subtree
          const auto min_d_travelled = calculate_d(start_pos, min_node);
          if(d < min_d_travelled) {
            // right node doesn't have sample. go back to left child
            node_pos -= 1;
            cost -= tree[node_pos].first;
            continue;
          }

          current_height += 1;
          cost += tree[node_pos << 1].first;
          node_pos = (node_pos << 1) + 1;
        }
        else {
          // if not good then go to left child
          current_height += 1;
          node_pos = node_pos << 1;
        }
      }
    }

    // check min_d_travelled before descending to the next node
    // if its < 0 then that means we have moved past the interval our sample is contained in.
    // if this is the case go to our left sibling
    ulint min_node = (node_pos << (height - current_height));
    const ulint min_d_travelled = calculate_d(start_pos, min_node);
    if(d < min_d_travelled) {
      node_pos -= 1;
      cost -= tree[node_pos].first;
      descend(start_pos, node_pos, cost, d, current_height); // descend again just in case there is tree left to be traversed
    }
  }

  //! Calculates the distance (number of leaf nodes between the two) between two node positions.
  /*!
    \param start_pos Left node to count from.
    \param end_pos Right node to count until.
  */
  ulint calculate_d(ulint start_pos, ulint end_pos) {
    ulint d = 0;
    d += leaf_node_bv.rank(end_pos);
    d -= leaf_node_bv.rank(start_pos);

    if(start_pos >= left_most_i) { // starting from bottom layer
      d -= leaf_node_bv.rank((start_pos >> 1));
      if(end_pos >= left_most_i) { // going to bottom layer bot->bot
        d += leaf_node_bv.rank((end_pos >> 1));
      }
      else { // going to top layer bot->top
        d += leaf_node_bv.rank((end_pos << 1));
      }
    } else { // starting from top layer
      d -= leaf_node_bv.rank((start_pos << 1));
      if(end_pos >= left_most_i) {// going to bottom layer top->bot
        d += leaf_node_bv.rank((end_pos >> 1));
      }
      else { // going to top layer top->top
        d += leaf_node_bv.rank((end_pos << 1));
      }
    }
    return d;
  }

  void print_tree_info() {
    cout << "tree size: " << tree.size() << endl;
    cout << "# of samples: " << leaf_samples.size() << endl;
    print_bounds();
    print_samples();
    print_nodebv();
    cout << "left_most_i: " << left_most_i << endl;
    cout << "height: " << height << endl;
  }

  void print_bounds() {
    cout << "tree bounds: ";
    for (size_t i = 0; i < tree.size(); i++) {
      cout << "[" << tree[i].first << "," << tree[i].second << "] ";
    }
    cout << endl;
  }

  void print_samples() {
    cout << "tree samples: ";
    for(size_t i = 0; i < leaf_samples.size(); i++) {
      cout << "[" << leaf_samples[i] << "] ";
    }
    cout << endl;
  }

  void print_nodebv() {
    cout << "node bv: ";
    for(size_t i = 0; i < leaf_node_bv.size(); i++) {
      cout << "[" << leaf_node_bv[i] << "] ";
    }
    cout << endl;
  }

  // serialize tree, leaf_samples, leaf_node_bv, left_most_i, height
  ulint serialize(std::ostream& out) {
    ulint w_bytes = 0;

    w_bytes += sdsl::serialize(tree.size(), out);
    out.write((char*)tree.data(), tree.size()*sizeof(tree[0]));
    w_bytes += sizeof(tree[0])*tree.size();

    w_bytes += sdsl::serialize(leaf_samples.size(), out);
    out.write((char*)leaf_samples.data(), leaf_samples.size()*sizeof(leaf_samples[0]));
    w_bytes += sizeof(leaf_samples[0])*leaf_samples.size();

    w_bytes += leaf_node_bv.serialize(out);
    w_bytes += sdsl::serialize(left_most_i, out);
    w_bytes += sdsl::serialize(height, out);

    return w_bytes;
  }

  // load tree, leaf_samples, leaf_node_bv, left_most_i, height
  void load(std::istream& in) {
    size_t temp_size;

    in.read((char*)&temp_size, sizeof(temp_size));
    tree.resize(temp_size);
    in.read((char*)tree.data(), temp_size*sizeof(tree[0]));

    in.read((char*)&temp_size, sizeof(temp_size));
    leaf_samples.resize(temp_size);
    in.read((char*)leaf_samples.data(), temp_size*sizeof(leaf_samples[0]));

    leaf_node_bv.load(in);

    in.read((char*)&left_most_i, sizeof(left_most_i));
    in.read((char*)&height, sizeof(height));
  }
};
}

#endif /* CSATREE_CLASSIC_HPP */
