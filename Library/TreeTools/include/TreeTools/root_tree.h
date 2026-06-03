#ifndef TreeTools_root_tree_
#define TreeTools_root_tree_

#include <Rcpp/Lightest>
#include <memory> /* for std::unique_ptr, make_unique */
#include <stdexcept> /* for errors */
#include "assert.h" /* for ASSERT */
#include "renumber_tree.h"
#include "types.h"

namespace TreeTools {
  extern inline Rcpp::IntegerMatrix preorder_edges_and_nodes(
      const Rcpp::IntegerVector parent,
      const Rcpp::IntegerVector child);

  extern inline std::pair<Rcpp::IntegerMatrix, Rcpp::NumericVector> preorder_weighted_pair(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::DoubleVector& weights);

  // edge must be BINARY
  // edge must be in preorder
  // 
  // Benchmarking at 2024-02-23 established that this is consistently twice
  // as fast as root_on_node, so is worth retaining,
  // despite some overlap in code.
  // 
  // #TODO Write test cases
  // 
  // [[Rcpp::export]]
  inline Rcpp::IntegerMatrix root_binary(const Rcpp::IntegerMatrix edge,
                                         const int outgroup) {

    const intx n_edge = edge.nrow(),
      n_node = n_edge / 2,
      n_tip = n_node + 1,
      root_node = n_tip + 1,
      max_node = n_node + n_tip;
    
    if (!n_edge || !n_node || n_tip < 2) return edge;
    if (edge(0, 1) == outgroup) return edge;
    if (outgroup < 1) {
      Rcpp::stop("`outgroup` must be a positive integer");
    }
    if (outgroup > max_node) {
      Rcpp::stop("`outgroup` exceeds number of nodes");
    }
    if (outgroup == root_node) {
      return edge;
    }
    

    std::unique_ptr<intx[]> edge_above = std::make_unique<intx[]>(max_node + 1);
    intx root_edges[2] = {0, 0};

    for (intx i = n_edge; i--; ) {

      edge_above[edge(i, 1)] = i;

      if (edge(i, 0) == root_node) {
        if (edge(i, 1) == outgroup) {
          return edge;
        }
        root_edges[root_edges[1] ? 0 : 1] = i;
      }

    }

    Rcpp::IntegerMatrix ret = Rcpp::clone(edge);
    intx invert_next = edge_above[outgroup];

    // We'll later add an edge from the now-unallocated root node to the outgroup.
    ret(invert_next, 0) = root_node;
    ret(invert_next, 1) = edge(invert_next, 0);

    do {
      invert_next = edge_above[edge(invert_next, 0)];
      ret(invert_next, 0) = edge(invert_next, 1);
      ret(invert_next, 1) = edge(invert_next, 0);
    } while (edge(invert_next, 0) != root_node);

    // second root i.e. 16 -- 24 must be replaced with root -> outgroup.
    intx spare_edge = (ret(root_edges[0], 0) == root_node ? 0 : 1);
    ret(invert_next, 1) = edge(root_edges[spare_edge], 1);
    ret(root_edges[spare_edge], 1) = outgroup;

    return preorder_edges_and_nodes(ret(Rcpp::_, 0), ret(Rcpp::_, 1));
  }

  // #TODO Write test cases
  // NB: If specifying internal node by number, note that node numbers will
  // change if tree is not already in preorder.
  // NB: root_node must == n_tip + 1
  //
  //  [[Rcpp::export]]
  inline Rcpp::List root_on_node(const Rcpp::List phy, const int outgroup) {

    Rcpp::IntegerMatrix edge = phy["edge"];
    Rcpp::NumericVector weight;

    const intx n_edge = edge.nrow();
    const intx n_node = phy["Nnode"];
    const intx max_node = n_edge + 1;
    const intx n_tip = max_node - n_node;
    const intx root_node = n_tip + 1;
    
    if (!n_edge || !n_node || n_tip < 2) return phy;

    const bool weighted = phy.containsElementNamed("edge.length");
    if (weighted) {
      
      Rcpp::IntegerVector parent_col(edge(Rcpp::_, 0));
      Rcpp::IntegerVector child_col(edge(Rcpp::_, 1));
      Rcpp::NumericVector edge_len = phy["edge.length"];
      std::tie(edge, weight) = preorder_weighted_pair(
        parent_col,
        child_col,
        edge_len
      );
    } else {
      edge = preorder_edges_and_nodes(edge(Rcpp::_, 0), edge(Rcpp::_, 1));
    }
    if (outgroup < 1) {
      Rcpp::stop("`outgroup` must be a positive integer");
    }
    if (outgroup > max_node) {
      Rcpp::stop("`outgroup` exceeds number of nodes");
    }
    Rcpp::List ret = Rcpp::clone(phy);
    if (outgroup == root_node) {
      ret.attr("order") = "preorder"; /* by preorder_weighted or _edges_&_nodes */
      ret["edge"] = edge;
      if (weighted) {
        ret["edge.length"] = weight;
      }
      return ret;
    }


    auto edge_above = std::make_unique<intx[]>(max_node + 1);
    intx root_edges[] = {0, 0};
    intx root_edges_found = 0;

    for (intx i = n_edge; i--; ) {
      edge_above[edge(i, 1)] = i;
      if (edge(i, 0) == root_node) {
        if (root_edges_found < 2) {
          root_edges[root_edges_found] = i;
        }
        ++root_edges_found;
      }
    }

    intx invert_next = edge_above[outgroup];
    
    
    if (root_edges_found == 2) { // Root node is vapour, and can be repurposed
      
      if (edge(root_edges[0], 1) == outgroup ||
          edge(root_edges[1], 1) == outgroup) {
        return phy;
      }
      
      Rcpp::IntegerVector new_parent = edge(Rcpp::_, 0);
      Rcpp::IntegerVector new_child = edge(Rcpp::_, 1);
      
      // We'll later add an edge from the now-unallocated root node to the outgroup.
      new_parent[invert_next] = root_node;
      new_child[invert_next] = edge(invert_next, 0);
      
      do {
        invert_next = edge_above[edge(invert_next, 0)];
        new_parent[invert_next] = edge(invert_next, 1);
        new_child[invert_next] = edge(invert_next, 0);
      } while (edge(invert_next, 0) != root_node);
      
      // Further root edges must be replaced with root -> outgroup
      intx spare_edge = (new_parent[root_edges[0]] == root_node ? 0 : 1);
      new_child[invert_next] = edge(root_edges[spare_edge], 1);
      new_child[root_edges[spare_edge]] = outgroup;
      
      if (weighted) {
        std::tie(edge, weight) = preorder_weighted_pair(new_parent, new_child, weight);
        ret["edge"] = edge;
        ret["edge.length"] = weight;
      } else {
        ret["edge"] = preorder_edges_and_nodes(new_parent, new_child);
      }
    } else { // Root node will be retained; we need a new root edge
      
      Rcpp::IntegerVector new_parent = edge(Rcpp::_, 0);
      Rcpp::IntegerVector new_child = edge(Rcpp::_, 1);
      
      const intx new_root = max_node + 1;
      
      new_parent.push_back(new_root);
      new_child.push_back(outgroup);

      new_parent[invert_next] = new_root;
      new_child[invert_next] = edge(invert_next, 0);

      while (edge(invert_next, 0) != root_node) {
        invert_next = edge_above[edge(invert_next, 0)];
        new_parent[invert_next] = edge(invert_next, 1);
        new_child[invert_next] = edge(invert_next, 0);
      }

      ret["Nnode"] = n_node + 1;
      if (weighted) {
        Rcpp::NumericVector new_wt(n_edge + 1);
        std::copy(weight.begin(), weight.end(), new_wt.begin());
        new_wt[n_edge] = 0;
        std::tie(edge, weight) = preorder_weighted_pair(new_parent, new_child, new_wt);
        ret["edge"] = edge;
        ret["edge.length"] = weight;
      } else {
        ret["edge"] = preorder_edges_and_nodes(new_parent, new_child);
      }
      
    }
    ret.attr("order") = "preorder"; /* by preorder_weighted or _edges_&_nodes */
    // #TODO there is probably a clever way to avoid doing a full preorder rewriting.
    return ret;
  }
}

#endif
