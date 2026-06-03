#ifndef TreeTools_n_cherries_
#define TreeTools_n_cherries_

#include <memory>    /* for std::unique_ptr */
#include <stdexcept> /* for errors */
#include <vector>

#include "assert.h" /* for ASSERT */

namespace TreeTools{

// Number of cherries in a binary phylogenetic tree
inline int n_cherries(const int* parent,
                      const int* child,
                      const size_t n_edge,
                      const int n_tip) {
  
  const size_t n_node = n_edge / 2;
  std::unique_ptr<bool[]> internal(new bool[n_node]());
  
  const bool unrooted = n_edge % 2;
  if (unrooted) {
    std::unique_ptr<bool[]> is_child(new bool[n_node + n_tip + 1]());

    for (size_t ed = 0; ed < n_edge; ++ed) {
      is_child[child[ed]] = true;
    }
    
    const int i_limit = n_tip + n_node + 1;
    int root_node = n_tip + 1;
    for (; root_node <= i_limit; ++root_node) {
      if (!is_child[root_node]) break;
    }
    
    if (root_node == i_limit) {
      throw std::runtime_error("Tree must be acyclic"); // nocov
    }
    
    bool root_internal_found = false;
    for (size_t ed = 0; ed < n_edge; ++ed) {
      const int child_i = child[ed];
      if (child_i > n_tip) {
        const int node = parent[ed];
        if (!root_internal_found && node == root_node) {
          root_internal_found = true;
        } else {
          internal[node - n_tip] = true;
        }
      }
    }
    
  } else {
    for (size_t ed = 0; ed < n_edge; ++ed) {
      if (child[ed] > n_tip) {
        const size_t node_idx = parent[ed] - n_tip;
        internal[node_idx] = true;
      }
    }
  }
  
  int n_cherries = 0;
  for (size_t i = 0; i < n_node; ++i) {
    if (!internal[i]) ++n_cherries;
  }
  return n_cherries;
}

}

#endif
