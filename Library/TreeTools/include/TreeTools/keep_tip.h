#ifndef TreeTools_keep_tip_
#define TreeTools_keep_tip_

#include <Rcpp/Lightest>
#include <stdexcept> /* for errors */
#include <memory> /* for make_unique */
#include <limits> /* for INT_MAX */
#include "assert.h" /* for ASSERT */
// #define TTDEBUG

#define RCOUT Rcpp::Rcout

#define RETAIN 9000

#define GET_NEW_NO(n) if (!new_no[n]) new_no[n] = ++next_no

#define SKIP_EDGE                                              \
GET_NEW_NO(parent);                                            \
new_no[child] = new_no[parent];                                \
if (rm_root && parent <= root_placeholder) {                   \
  root_placeholder = child;                                    \
}


namespace TreeTools {
  // entry 0, 1... of keep is TRUE if leaf "1", "2"... should be retained, 
  // false otherwise.
  // [[Rcpp::export]]
  inline Rcpp::IntegerMatrix keep_tip (
      const Rcpp::IntegerMatrix edge,
      const Rcpp::LogicalVector keep) {
    
    if (edge.ncol() != 2) {
      Rcpp::stop("edge must have two columns");
    }
    
    ASSERT(keep.length() < INT_MAX - 1);
    const int
      start_edge = edge.nrow(),
      n_tip = int(keep.length()),
      root_node = n_tip + 1,
      all_nodes = start_edge + 2
    ;
    
    auto
      n_child = std::make_unique<int[]>(all_nodes),
      one_child = std::make_unique<int[]>(all_nodes),
      new_no = std::make_unique<int[]>(all_nodes)
    ;
    
    int next_no = 0;
    for (int i = 0; i != n_tip; ++i) {
      if (keep[i]) {
        n_child[i + 1] = RETAIN;
  #ifdef TTDEBUG
        n_child[i + 1] = RETAIN + i + 1;
  #endif
        new_no[i + 1] = ++next_no;
      }
    }
    
    // *** Initialize postorder traversal *** //
    int
      kept_edges = 0,
        root_edges = 0
    ;
    for (int i = start_edge; i--; ) {
      const int
      edge_parent = edge(i, 0),
        edge_child = edge(i, 1),
        downstream = n_child[edge_child]
      ;
      if (edge_parent == root_node) {
        ++root_edges;
      }
      if (downstream) {
        ++n_child[edge_parent];
        if (downstream == 1) {
          one_child[edge_parent] = one_child[edge_child];
        } else {
          one_child[edge_parent] = edge_child;
          ++kept_edges;
        }
      }
    }
    int new_root = root_node;
    if (n_child[root_node] == 1) {
      // We've deleted one side of the tree entirely, creating a trailing root
      --kept_edges;
      new_root = one_child[root_node];
    }
    
  #ifdef TTDEBUG
    RCOUT << " > New basal node: " << new_root << ".\n";
    for (int i = 0; i != start_edge; ++i) {
      RCOUT << " - Edge " << (1+i) << ": " << n_child[edge(i, 1)] << ".\n";
    }
  #endif
    
    // Keep unrooted trees unrooted
    const bool rooted = root_edges == 2;
    bool rm_root = !rooted && n_child[new_root] == 2;
  #ifdef TTDEBUG
    RCOUT << " > Originally " << (rooted ? "" : "un") << "rooted; "
          << "new root " << new_root << " has " << n_child[new_root]
          << " surviving children.\n\n";
  #endif
    if (rm_root) {
      --kept_edges;
    }
    
    // *** Initialize preorder traversal *** //
    int writing_edge = -1;
    int root_placeholder = root_node;
    Rcpp::IntegerMatrix ret(kept_edges, 2);
    
    for (int i = 0; i != start_edge; ++i) {
      const int
      parent = edge(i, 0),
        child = edge(i, 1),
        n_children = n_child[child]
      ;
      if (n_children) {
        if (n_children == 1) {
          SKIP_EDGE;
  #ifdef TTDEBUG
          RCOUT << " x Omit " << edge(i, 0) << "-" << edge(i, 1)
                << ": \"" << new_no[parent] << "\"\n";
  #endif
        } else {
          if (n_child[root_node] == 1) {
            // Skip this dangling root edge, but mark root as ready to retain
  #ifdef TTDEBUG
            RCOUT << " x Skip dangling root edge " << edge(i, 0) << "-"
                  << edge(i, 1) << ".\n";
  #endif
            SKIP_EDGE;
            n_child[root_node] = RETAIN;
          } else if (rm_root
                       && parent <= root_placeholder
                       && child > n_tip) {
                       // This is the duplicated root edge
                       rm_root = false; // We will not write it
            GET_NEW_NO(parent);
            new_no[child] = new_no[parent];
  #ifdef TTDEBUG
            RCOUT << " x Skip duplicated root edge " << edge(i, 0) << "-"
                  << edge(i, 1) << " (\"" << new_no[child] << "\").\n";
  #endif
          } else {
            // Record this edge:
            ++writing_edge;
            ASSERT(writing_edge < kept_edges);
            GET_NEW_NO(parent);
            ret(writing_edge, 0) = new_no[parent];
            GET_NEW_NO(child);
            ret(writing_edge, 1) = new_no[child];
#ifdef TTDEBUG
            RCOUT << " > Translating: " 
                  << edge(i, 0) << "-" << edge(i, 1)
                  << " --> " << new_no[parent] << "-" << new_no[child] << "\n";
#endif
          }
        }
      }
    }
    
    return ret;
  }
}

#endif
