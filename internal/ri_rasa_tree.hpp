#ifndef RI_RASA_TREE_HPP
#define RI_RASA_TREE_HPP

#include "definitions.hpp"
#include "r_index.hpp"
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"
#include "csatree_classic.hpp"
#include "dcheck.hpp"
#include "vector.hpp"

using namespace sdsl;
// #define BISIMULATE_CSA_TREE 1
#define WITH_LEAF_NODE_BV 1

namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>

//! Tree data structure that stores the samples contained within a cycle. Performs queries.
class rads_tree {

#ifdef BISIMULATE_CSA_TREE
  CSATree<sparse_bv_type, rle_string_t> m_classic_tree; //DEBUG
  std::vector<std::tuple<ulint, ulint, uint>> m_tree_pointers;//DEBUG
#endif//BISIMULATE_CSA_TREE_CSA_TREE

public:
  static constexpr size_t VERSION = 1;
  using cost_type = ulint;
  using limit_type = long long int;

  Vector<std::pair<cost_type, limit_type>> m_tree; // nodes are pairs that represent: (edge cost, edge threshold).
  Vector<ulint> leaf_samples; //@ stores the SA sample values of the leaves represented by run indices
#ifdef WITH_LEAF_NODE_BV 
  sparse_bv_type leaf_node_bv; // bv telling us which node is a leaf node //TODO: can be removed!
#endif//WITH_LEAF_NODE_BV
  uint left_most_i; // index at which the left most leaf is stored || sdsl::serialize
  size_t height; // height of the tree starting at 0 || sdsl::serialize

  rads_tree(){};
  //! rads_tree constructor that recursively builds a balanced binary interval tree using the samples provided.
  /*!
    \param cycle Path of sample nodes on which to be built a tree.
    \param bounds Calculated bounds used for lower and upper bound calculation.
    \param tree_num Number of tree.
    \param tree_pointers Vector of tree pointers that is added to every time a leaf is set.
  */
  rads_tree(const std::vector<ulint> &cycle, const std::vector<std::pair<ulint, ulint>> &bounds, const uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers) {
    size_t n = cycle.size(); // path size

    height = ceil(log2(n)); // height of the tree //TODO: can be removed
    //TODO: make std::move(cycle)

#ifdef BISIMULATE_CSA_TREE
    m_classic_tree = decltype(m_classic_tree)(cycle, bounds, tree_num, m_tree_pointers);
#endif//BISIMULATE_CSA_TREE_CSA_TREE
    leaf_samples = cycle; // store the samples in the cycle


    const size_t number_of_perfect_leaves = 1ULL<<height;
    DCHECK_LE(n, number_of_perfect_leaves);
    const size_t round_off = number_of_perfect_leaves - n; //@ number of excess nodes from the next power of two
    const size_t number_of_perfect_nodes = (1ULL<<(height + 1));
    m_tree = decltype(m_tree)(number_of_perfect_nodes - round_off, std::make_pair(0, -1)); // tree size initialization, default the costs to -1

    auto temp_leaf_bv = vector<bool>(number_of_perfect_nodes, false); //TODO

    constructor_helper(bounds, 1, 0, number_of_perfect_leaves - 1, tree_num, tree_pointers, temp_leaf_bv);

    left_most_i = (number_of_perfect_nodes>>1);

#ifdef WITH_LEAF_NODE_BV 
    leaf_node_bv = sparse_bv_type(temp_leaf_bv);
    for(size_t i = 0; i < left_most_i; ++i) {
      DCHECK_EQ(leaf_node_bv[i], false);
    }
    for(size_t i = left_most_i; i < m_tree.size(); ++i) {
      DCHECK_EQ(leaf_node_bv[i], true);
    }
#endif//WITH_LEAF_NODE_BV
  }

  rads_tree(const rads_tree &other_tree) {
    this->m_tree = other_tree.m_tree;
    this->leaf_samples = other_tree.leaf_samples;
#ifdef WITH_LEAF_NODE_BV 
    this->leaf_node_bv = other_tree.leaf_node_bv;
#endif//WITH_LEAF_NODE_BV
    this->left_most_i = other_tree.left_most_i;
    this->height = other_tree.height;
#ifdef BISIMULATE_CSA_TREE
    this->m_classic_tree = other_tree.m_classic_tree;
    this->m_tree_pointers = other_tree.m_tree_pointers;
#endif//BISIMULATE_CSA_TREE_CSA_TREE
  }

  rads_tree(rads_tree &&other_tree) noexcept
  : m_tree(move(other_tree.m_tree))
  , leaf_samples(move(other_tree.leaf_samples))
#ifdef WITH_LEAF_NODE_BV 
  , leaf_node_bv(move(other_tree.leaf_node_bv))
#endif//WITH_LEAF_NODE_BV
  , left_most_i(move(other_tree.left_most_i))
  , height(move(other_tree.height))
#ifdef BISIMULATE_CSA_TREE
  , m_classic_tree(move(other_tree.m_classic_tree))
  , m_tree_pointers(move(other_tree.m_tree_pointers))
#endif//BISIMULATE_CSA_TREE_CSA_TREE
  {}

  rads_tree& operator=(const rads_tree &other_tree) {
    return *this = rads_tree(other_tree);
  }

  rads_tree& operator=(rads_tree &&other_tree) {
    swap(m_tree, other_tree.m_tree);
    swap(leaf_samples, other_tree.leaf_samples);
#ifdef WITH_LEAF_NODE_BV 
    swap(leaf_node_bv, other_tree.leaf_node_bv);
#endif//WITH_LEAF_NODE_BV
    swap(left_most_i, other_tree.left_most_i);
    swap(height, other_tree.height);
#ifdef BISIMULATE_CSA_TREE
    swap(m_classic_tree, other_tree.m_classic_tree);
    swap(m_tree_pointers, other_tree.m_tree_pointers);
#endif//BISIMULATE_CSA_TREE_CSA_TREE
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
    if(begin >= end) { // condition which creates the leaf nodes.
      if(node < m_tree.size()) {
        m_tree[node].first = bounds[leaf_samples[begin]].first; // edge cost
        m_tree[node].second = bounds[leaf_samples[begin]].second; // edge threshold
        if(begin != leaf_samples.size()-1) { // if it is not the last sample in the tree, add a tree pointer to be used when querying
          tree_pointers.push_back(std::make_tuple(leaf_samples[begin], tree_num, node));
        }
        leaf_bv[node] = true;
      }
      return;
    }
    const size_t mid = (begin + end)/2;
    DCHECK_EQ(mid-begin, end-mid-1);

    constructor_helper(bounds, 2*node, begin, mid, tree_num, tree_pointers, leaf_bv); // construct left child
    constructor_helper(bounds, 2*node+1, mid+1, end, tree_num, tree_pointers, leaf_bv); // construct right child
    if(2*node+1 < m_tree.size()) {
      m_tree[node].first = m_tree[2*node].first + m_tree[2*node+1].first; // edge sum of left and right child
      m_tree[node].second = std::min<limit_type>(m_tree[2*node].second, m_tree[2*node+1].second - m_tree[2*node].first); // edge threshold
    }
  }


bool can_traverse(const ulint node_id, const cost_type cost) const {
    return node_id < m_tree.size() 
      && m_tree[node_id].second >= 0  //@ costs are signed, so first check whether we can scale it up to unsigned
      && m_tree[node_id].second > cost;
}

ulint lowest_left_ancestor(const ulint node_id) const {
      return node_id >> (__builtin_ctzl(~node_id));
}

ulint node_height(const ulint node_id) const {
  DCHECK_GT(node_id, 0);
  DCHECK_LT(sdsl::bits::hi(node_id),  tree_height());
  return tree_height() - sdsl::bits::hi(node_id)-1;
}
uint tree_height() const {
  DCHECK_GT(m_tree.size(), 0);
    return sdsl::bits::hi(m_tree.size()-1)+1;
}

ulint leftmost_leaf() const {
  const auto ret = 1ULL << (tree_height()-1);
  DCHECK_EQ(ret, left_most_i);
  return ret;
}

ulint leftmost_leaf(const ulint node_id) const {
  return node_id<<node_height(node_id);
}

bool is_leaf(const ulint node_id) const {
  return node_id >= leftmost_leaf();
}

ulint number_of_leaves(const ulint node_id) const {
  return 1ULL<<node_height(node_id);
}


std::tuple<ulint, ulint, ulint>  query(const ulint leaf_pos, uint cost, uint distance_bound, ulint = -1) const {
  if(leaf_pos + distance_bound >= m_tree.size()) { distance_bound = m_tree.size()-leaf_pos-1;  } //@ no overshooting
  uint initial_cost = cost;
  ulint visiting_node = leaf_pos;

  auto remaining_distance = [&] () {
    return distance_bound - (leftmost_leaf(visiting_node) - leaf_pos);
  };
  
  if(!can_traverse(visiting_node, cost)) {
    DCHECK(false);
  }
  // visiting_node+=1;
  cost += m_tree[leaf_pos].first; 
  if( (visiting_node & 1) == 0) {
    visiting_node += 1;
  } //! invariant: always start from a right child being a leaf

  bool ascended = false;
  //while(can_traverse(visiting_node, cost)) {
  while(true) {
    const auto left_ancestor = lowest_left_ancestor(visiting_node);
    if(left_ancestor == 0) {
      visiting_node = 1;
      break;
    }
    visiting_node = left_ancestor;
    if(!can_traverse(visiting_node, cost) || number_of_leaves(visiting_node) > remaining_distance()) {
      break;
    }
    cost += m_tree[visiting_node+1].first; // take the costs of our right sibling
    visiting_node >>= 1; //parent node
                         ascended = true;
  }
  DCHECK_NE(visiting_node, leaf_pos);
  DCHECK_LT(node_height(leaf_pos), node_height(visiting_node));
  // if we cannot further climb up, go down again and to the right child
  visiting_node <<=1;
  visiting_node += 1;
  if(ascended) {
  cost -= m_tree[visiting_node].first; // revert adding the cost of our last right sibling
  }

  DCHECK_NE(lowest_left_ancestor(visiting_node), 1); //cannot go to right
  while(can_traverse(visiting_node, cost) && node_height(visiting_node) == node_height(visiting_node+1) && number_of_leaves(visiting_node) <= remaining_distance()) {
    cost += m_tree[visiting_node].first;
    visiting_node += 1;
  }

  while(!is_leaf(visiting_node)) {
    visiting_node <<= 1; //left child node 
    DCHECK_EQ(visiting_node & 1, 0);
    if(can_traverse(visiting_node, cost) && number_of_leaves(visiting_node) <= remaining_distance()) { // go to right child if possible
      cost += m_tree[visiting_node].first;
      visiting_node += 1; // right sibling
      DCHECK_EQ(visiting_node & 1, 1);
    }
  }
  ulint distance = visiting_node - leaf_pos;
  const ulint leftmost_leaf_pos = this->left_most_i;

  if(distance > distance_bound) { distance = distance_bound; }
  if(leaf_pos + distance >= m_tree.size()) { distance = m_tree.size()-1-leaf_pos; }

  const auto ret = std::make_tuple(this->leaf_samples[distance+(leaf_pos - leftmost_leaf_pos)], distance, cost); 

  cost -= m_tree[visiting_node].first; //@ do not charge for the last node

  {
  const ulint leaf_distance = leaf_query(leaf_pos, initial_cost, distance_bound);
  const auto myret = std::make_tuple(this->leaf_samples[leaf_distance+(leaf_pos - leftmost_leaf_pos)], leaf_distance, initial_cost); 
  check_query_tuple(ret, myret);
  }
  return ret;
}

static void check_query_tuple(const std::tuple<ulint, ulint, ulint>& a, const std::tuple<ulint, ulint, ulint>& b)  {
  DCHECK_EQ(std::get<0>(a), std::get<0>(b));
  DCHECK_EQ(std::get<1>(a), std::get<1>(b));
  DCHECK_EQ(std::get<2>(a), std::get<2>(b));
}

ulint leaf_query(const ulint leaf_pos, uint& cost, const uint distance_bound) const {
  const auto& treearray =  this->m_tree;
  const ulint leftmost_leaf_pos = this->left_most_i;
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
  return distance;
}


  //! Query function.
  /*!
    \param node_pos Leaf node to begin traversal.
    \param cost Cost that we will carry and add to as we climb.
    \param d Distance that we would ideally travel.
  */
  std::tuple<ulint, ulint, ulint> tree_query(ulint node_pos, uint cost, uint d, 
#ifdef BISIMULATE_CSA_TREE
      const ulint leaf_sample = -1
#else
      const ulint = -1
#endif
      ) {
    //! Climbing traversal function. Climb until we cannot and then begin descending.
    // begin climbing the tree and only collect when we go to a sibling.
    // when you go to a parent do not collect the cost, that is if the parent has the same left most node.
    // if youre moving to a subtree where the leftmost node is not the same then we have to collect the money.
   

// #ifdef BISIMULATE_CSA_TREE
//     const auto initial_cost = cost;
//     const auto initial_d = d;
// #endif//BISIMULATE_CSA_TREE_CSA_TREE
    const ulint start_pos = node_pos;

    ulint current_height = height;
    ulint prev_pos = node_pos;
    ulint max_d_travelled = 0;

    // we can enter at either the height or height - 1
    if(node_pos < left_most_i) {
      DCHECK(false); // cannot happen!
      // current_height -= 1;
    }

    // climb up while the upper_bounds let us / climb up at least once if we are still at the start node
    ulint current_leftmost_node = (node_pos << (height - current_height));
    while(((node_pos != 1) && (m_tree[node_pos].second > 0) && (cost < m_tree[node_pos].second)) || (start_pos == node_pos)) {
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
        cost += m_tree[prev_pos].first;
        descend(start_pos, node_pos, cost, d, current_height);
        if(node_pos >= m_tree.size()) {
          node_pos--;
          cost -= m_tree[node_pos].first;
          DCHECK_LT(node_pos, m_tree.size());
        }
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
          if((m_tree[node_pos >> 1].second > 0) && (cost < m_tree[node_pos >> 1].second)) {
            node_pos = (node_pos >> 1) + 1;
          }
          else {
            break; // get out cause you're going to make an illegal jump
          }
        }
        else {
          node_pos = (node_pos >> 1) + 1;
          DCHECK_LT(node_pos, m_tree.size());
        }
      }

      current_height -= 1;

      // moved to a subtree where the left most node is not the same, therefore we collect
      if((node_pos << (height - current_height)) != current_leftmost_node) {
        cost += m_tree[prev_pos].first;
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
      while(node_pos < left_most_i) { //@ while node is an internal node
      // while(node_pos < leaf_node_bv.size() && !leaf_node_bv[node_pos]) {
        if(node_pos >= m_tree.size() || !(m_tree[node_pos].second > 0) || !(cost < m_tree[node_pos].second)) { // this means we keep getting denied
          node_pos = (node_pos << 1);
          current_height += 1;
        }
        else {
          node_pos = (node_pos << 1) + 1; // once we find a bound that were allowed to go through, go to the right child.
          current_height += 1;
          cost += m_tree[node_pos - 1].first;
          descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
          break;
        }
      }
    }
    else {
      cost += m_tree[node_pos << 1].first; // add cost of node that denies us
      node_pos = std::min((node_pos << 1) + 1, m_tree.size()-1); // move node_pos to the right child of the current node // or should it be the left child of the node
      current_height += 1;
      descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
    }

    bool overshoot = false;
    if(node_pos >= m_tree.size()) {  // no overshooting
      node_pos = m_tree.size()-1; 
      cost -= m_tree[node_pos].first;
      overshoot = true;
    }

    ulint distance = calculate_d(start_pos, node_pos); // calculate the distance from start to end to know how far we travelled
    ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
    const auto my_ret = std::make_tuple(leaf_samples[leaf_sample_index], distance, cost); // return new sample and new distance


#ifdef BISIMULATE_CSA_TREE
    // if(leaf_sample != -1ULL) {
    //   bool found = false;
    //   for(size_t i = 0; i < m_tree_pointers.size(); ++i) {
    //     if(std::get<0>(m_tree_pointers[i]) == leaf_sample) {
    //       const auto ret = m_classic_tree.query(std::get<2>(m_tree_pointers[i]), initial_cost, initial_d);
    //       CHECK_EQ(std::get<0>(ret), std::get<0>(my_ret));
    //       CHECK_EQ(std::get<1>(ret), std::get<1>(my_ret));
    //       CHECK_EQ(std::get<2>(ret), std::get<2>(my_ret));
    //       found = true;
    //       break;
    //     }
    //   }
    //   CHECK_EQ(found, true);
    // }
#endif//BISIMULATE_CSA_TREE_CSA_TREE
    return my_ret;
  }

  //! Descent function. This descends until we reach a leaf node making the proper movements based on our cost and the bounds we encounter.
  /*!
    \param start_pos Node we started the whole traversal from.
    \param node_pos Node position that we are currently at.
    \param cost Cost we have accumulated so far.
    \param current_height Current height in the tree that we're at.
  */
  void descend(ulint start_pos, ulint &node_pos, uint &cost, uint d, ulint current_height) {
    DCHECK_LE(current_height, height);
    ulint prev_pos = node_pos;

    while(node_pos < m_tree.size()) { // while we are in the bounds of the tree and not at a leaf node
      const auto last_bit = (node_pos & 1);
      prev_pos = node_pos;
      if(last_bit == 0) { // left node descent
        if(node_pos < m_tree.size() && (m_tree[node_pos].second > 0) && (cost < m_tree[node_pos].second)) { // if good, go to right sibling
          const ulint min_node = ((node_pos + 1) << (height - current_height));
          const auto min_d_travelled = calculate_d(start_pos, min_node);
          if(d < min_d_travelled) {
            // this means we would go too far when going to right sibling
            // go to right child instead
            if(current_height == height) { break; }
            node_pos = node_pos << 1;
            node_pos += 1;
            DCHECK_LT(node_pos, m_tree.size());
            current_height += 1;
            cost += m_tree[node_pos - 1].first;
          }
          else {
            // go to right sibling
            if(node_pos + 1 == m_tree.size()) { break; } // already at the last leaf
            node_pos += 1;
            cost += m_tree[prev_pos].first;
          }
        }
        else { // if not good then go to left child and check on next iteration
          if(current_height == height) { break; }
          node_pos = node_pos << 1;
          // DCHECK_LT(node_pos, tree.size());
          current_height += 1;
        }
      }
      else { // right node descend
        if((m_tree[node_pos].second > 0) && (cost < m_tree[node_pos].second)) { // if good, go to right child
          const ulint min_node = (node_pos << (height - current_height)); // this is how you can get the left mode node of the current subtree
          const auto min_d_travelled = calculate_d(start_pos, min_node);
          if(d < min_d_travelled) {
            // right node doesn't have sample. go back to left child
            node_pos -= 1;
            cost -= m_tree[node_pos].first;
            continue;
          }

          if(current_height == height) { break; }
          current_height += 1;
          cost += m_tree[node_pos << 1].first;
          node_pos = (node_pos << 1) + 1;
          DCHECK_LT(node_pos, m_tree.size());
        }
        else {
          // if not good then go to left child
          if(current_height == height) { break; }
          current_height += 1;
          node_pos = node_pos << 1;
          // if(node_pos < m_tree.size()) { --node_pos; } // do not overshoot
          // DCHECK_LT(node_pos, m_tree.size());
        }
      }
    }

    // check min_d_travelled before descending to the next node
    // if its < 0 then that means we have moved past the interval our sample is contained in.
    // if this is the case go to our left sibling
    const ulint min_node = (node_pos << (height - current_height));
    DCHECK_LE(start_pos, min_node);
    const ulint min_d_travelled = calculate_d(start_pos, min_node);
    if(d < min_d_travelled) {
      node_pos -= 1;
      cost -= m_tree[node_pos].first;
      descend(start_pos, node_pos, cost, d, current_height); // descend again just in case there is tree left to be traversed
    }
  }

#ifdef WITH_LEAF_NODE_BV 
  ulint calculate_d_debug(ulint start_pos, ulint end_pos) {
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
#endif//WITH_LEAF_NODE_BV

  //! Calculates the distance (number of leaf nodes between the two) between two node positions.
  /*!
    \param start_pos Left node to count from.
    \param end_pos Right node to count until.
  */
  ulint calculate_d(ulint start_pos, ulint end_pos) {
    DCHECK_LE(start_pos, end_pos);
    const auto ret = std::min(end_pos,m_tree.size())-start_pos;

#ifdef WITH_LEAF_NODE_BV
    DCHECK_EQ(ret, calculate_d_debug(start_pos, end_pos)); //TODO
#endif//WITH_LEAF_NODE_BV
    return std::min(end_pos,m_tree.size())-start_pos;
  }

  void print_tree_info() {
    cout << "m_tree size: " << m_tree.size() << endl;
    cout << "# of samples: " << leaf_samples.size() << endl;
    print_bounds();
    print_samples();
#ifdef WITH_LEAF_NODE_BV 
    print_nodebv();
#endif//WITH_LEAF_NODE_BV
    cout << "left_most_i: " << left_most_i << endl;
    cout << "height: " << height << endl;
  }

  void print_bounds() {
    cout << "m_tree bounds: ";
    for (size_t i = 0; i < m_tree.size(); i++) {
      cout << "[" << m_tree[i].first << "," << m_tree[i].second << "] ";
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

#ifdef WITH_LEAF_NODE_BV 
  void print_nodebv() {
    cout << "node bv: ";
    for(size_t i = 0; i < leaf_node_bv.size(); i++) {
      cout << "[" << leaf_node_bv[i] << "] ";
    }
    cout << endl;
  }
#endif//WITH_LEAF_NODE_BV

  // serialize tree, leaf_samples, leaf_node_bv, left_most_i, height
  ulint serialize(std::ostream& out) {
    ulint w_bytes = 0;

    w_bytes += sdsl::serialize(m_tree.size(), out);
    out.write((char*)m_tree.data(), m_tree.size()*sizeof(m_tree[0]));
    w_bytes += sizeof(m_tree[0])*m_tree.size();

    w_bytes += sdsl::serialize(leaf_samples.size(), out);
    out.write((char*)leaf_samples.data(), leaf_samples.size()*sizeof(leaf_samples[0]));
    w_bytes += sizeof(leaf_samples[0])*leaf_samples.size();

#ifdef WITH_LEAF_NODE_BV 
    w_bytes += leaf_node_bv.serialize(out);
#endif//WITH_LEAF_NODE_BV
    w_bytes += sdsl::serialize(left_most_i, out);
    w_bytes += sdsl::serialize(height, out);

#ifdef BISIMULATE_CSA_TREE
    w_bytes += m_classic_tree.serialize(out);

    w_bytes += sdsl::serialize(m_tree_pointers.size(), out);
    out.write((char*)m_tree_pointers.data(), m_tree_pointers.size()*sizeof(m_tree_pointers[0]));
    w_bytes += sizeof(m_tree_pointers[0])*m_tree_pointers.size();
#endif//BISIMULATE_CSA_TREE_CSA_TREE
    return w_bytes;
  }

  // load tree, leaf_samples, leaf_node_bv, left_most_i, height
  void load(std::istream& in) {
    size_t temp_size;

    in.read((char*)&temp_size, sizeof(temp_size));
    m_tree.resize(temp_size);
    in.read((char*)m_tree.data(), temp_size*sizeof(m_tree[0]));

    in.read((char*)&temp_size, sizeof(temp_size));
    leaf_samples.resize(temp_size);
    in.read((char*)leaf_samples.data(), temp_size*sizeof(leaf_samples[0]));

#ifdef WITH_LEAF_NODE_BV 
    leaf_node_bv.load(in);
#endif//WITH_LEAF_NODE_BV

    in.read((char*)&left_most_i, sizeof(left_most_i));
    in.read((char*)&height, sizeof(height));

#ifdef BISIMULATE_CSA_TREE
    m_classic_tree.load(in);

    in.read((char*)&temp_size, sizeof(temp_size));
    m_tree_pointers.resize(temp_size);
    in.read((char*)m_tree_pointers.data(), temp_size*sizeof(m_tree_pointers[0]));
#endif//BISIMULATE_CSA_TREE_CSA_TREE
  }
};
}

#endif /* RI_RASA_TREE_HPP */
