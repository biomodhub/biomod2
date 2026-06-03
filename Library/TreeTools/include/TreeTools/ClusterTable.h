#ifndef _TREETOOLS_CLUSTERTABLE_H
#define _TREETOOLS_CLUSTERTABLE_H

#include <array> /* for array */
#include <bitset> /* for bitset */
#include <vector> /* for vector */
#include <sstream> /* for ostringstream */
#include <Rcpp/Lightest>
#include "assert.h" /* for ASSERT */
#include "types.h" /* for int32 */
#include "root_tree.h" /* for root_on_node */

#define UNINIT -999
#define INF TreeTools::INTX_MAX

#define CT_IS_LEAF(a) (a) <= n_tip

namespace TreeTools {

  // Use stack allocation for trees up to this size for optimal performance
  inline constexpr int_fast32_t ct_stack_threshold = 8192;
  // New increased limit with heap allocation
  inline constexpr int_fast32_t ct_max_leaves_heap = 100000;
  
  template <typename T>
  inline void resize_uninitialized(std::vector<T>& v, std::size_t n) {
    static_assert(std::is_trivial<T>::value, "Requires trivial type");
    if (n > v.size()) {
      v.reserve(n);
      v.insert(v.end(), n - v.size(), T{});
    } else {
      v.resize(n);
    }
  }
  
  // Dynamic bitset that uses stack allocation for small sizes, heap for large
  class DynamicBitset {
  private:
    std::bitset<ct_stack_threshold + 1> stack_bits;
    std::vector<unsigned char> heap_bits;
    bool use_heap;
    std::size_t size_;
    
  public:
    DynamicBitset() : use_heap(false), size_(0) {}
    
    void resize(std::size_t n) {
      size_ = n;
      if (n > ct_stack_threshold) {
        use_heap = true;
        heap_bits.assign(n, 0);
      } else {
        use_heap = false;
        stack_bits.reset();
      }
    }
    
    bool operator[](std::size_t idx) const {
      return use_heap ? heap_bits[idx] != 0 : stack_bits[idx];
    }
    
    void set(std::size_t idx, bool value = true) {
      if (use_heap) {
        heap_bits[idx] = value ? 1 : 0;
      } else {
        stack_bits[idx] = value;
      }
    }
    
    void reset() {
      if (use_heap) {
        std::fill(heap_bits.begin(), heap_bits.end(), 0);
      } else {
        stack_bits.reset();
      }
    }
    
    std::size_t size() const { return size_; }
  };
  
  class ClusterTable {
    
    struct ClusterRow {
      int32 L;
      int32 R;
    };
    
    int32 n_edge;
    int32 n_internal;
    int32 n_leaves;
    int32 n_shared = 0;
    int32 enumeration = 0;
    int32 v_j;
    int32 Tlen;
    int32 Tlen_short;
    int32 X_ROWS;
    std::vector<int32> internal_label;
    int32 *internal_label_ptr = nullptr;
    std::vector<int32> leftmost_leaf;
    std::vector<int32> T;
    std::size_t T_idx = 0;
    std::vector<int32> visited_nth;
    std::vector<ClusterRow> x_rows;
    // Dynamic bitset that uses stack allocation for small trees,
    // heap allocation for large trees
    DynamicBitset Xswitch;
    // Track number of set switches (excluding index 0)
    std::size_t xswitch_set_count = 0;
    

  public:
    ClusterTable(Rcpp::List); // i.e. PREPARE(T)

    [[nodiscard]] inline bool is_leaf(const int32 v) noexcept {
      return v <= n_leaves;
    }
    
    [[nodiscard]] inline const int32 edges() noexcept {
      return n_edge;
    }

    [[nodiscard]] inline const int32 leaves() noexcept {
      return n_leaves;
    }

    inline void ENTER(int32 v, int32 w) noexcept {
      T[T_idx++] = v;
      T[T_idx++] = w;
    }

    [[nodiscard]] inline int32 N() noexcept {
      return n_leaves;
    }

    [[nodiscard]] inline int32 M() noexcept {
      return n_internal;
    }

    inline void TRESET() noexcept {
      // This procedure prepares T for an enumeration of its entries,
      // beginning with the first entry.
      T_idx = 0;
    }

    inline void READT(int32 *v, int32 *w) {
      *v = T[T_idx++];
      *w = T[T_idx++];
    }

    inline void NVERTEX(int32 *v, int32 *w) noexcept {
      if (T_idx != static_cast<size_t>(Tlen)) {
        READT(v, w);
        v_j = *v;
      } else {
        *v = 0;
        *w = 0;
      }
    }

    inline void NVERTEX_short(int32 *v, int32 *w) noexcept {
      // Don't count all-tips or all-ingroup: vertices 0, ROOT, Ingp.
      if (T_idx != static_cast<size_t>(Tlen_short)) {
        READT(v, w);
        // v_j = *v; // Unneeded unless we go on to call LEFTLEAF
      } else {
        *v = 0;
        *w = 0;
      }
    }

    inline int32 LEFTLEAF() noexcept {
      // If NVERTEX has returned entry <vj, wj> in T, the leftmost leaf in the
      // subtree rooted at vj has entry <vk, wk> where k = j - wj.
      // This function procedure returns Vk as its value.
      return leftmost_leaf[v_j - 1];
    }

    inline void SET_LEFTMOST(int32 index, int32 val) noexcept {
      leftmost_leaf[index - 1] = val;
    }

    [[nodiscard]] inline int32 GET_LEFTMOST(int32 index) noexcept {
      return leftmost_leaf[index - 1];
    }

    // Procedures to manipulate cluster tables, per Table 4 of Day 1985.

    inline int32 ENCODE(const int v) noexcept {
      // This function procedure returns as its value the internal label
      // assigned to leaf v
      // MS note: input = v; output = X[v, 3]
      return internal_label[v];
    }

    inline int32 DECODE(const int32 internal_relabeling) noexcept {
      // MS: input = X[v, 3], output = v
      return visited_nth[internal_relabeling - 1];
    }

    inline void VISIT_LEAF (const int32* leaf, int32* n_visited) {
      visited_nth[(*n_visited)++] = *leaf;
      internal_label[*leaf] = *n_visited;
    }

    Rcpp::IntegerVector X_decode() {
      Rcpp::IntegerVector ret(N());
      for (int32 i = n_leaves; i--; ) {
        ret(i) = DECODE(i + 1);
      }
      return ret;
    }
    
    [[nodiscard]] inline int32 X_left(int32 row) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      ASSERT(x_rows[row - 1].L < std::numeric_limits<int32>::max());
      return x_rows[row - 1].L;
    }
    
    [[nodiscard]] inline int32 X_right(int32 row) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      ASSERT(x_rows[row - 1].R < std::numeric_limits<int32>::max());
      return x_rows[row - 1].R;
    }
    
    inline void setX_left(int32 row, int32 value) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      x_rows[row - 1].L = value;
    }
    
    inline void setX_right(int32 row, int32 value) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      x_rows[row - 1].R = value;
    }
    
    Rcpp::IntegerMatrix X_contents() noexcept {
      Rcpp::IntegerMatrix ret(X_ROWS, 2);
      for (int32 i = X_ROWS; i--; ) {
        ret(i, 0) = x_rows[i].L;
        ret(i, 1) = x_rows[i].R;
      }
      return ret;
    }

    [[nodiscard]] inline bool CLUSTONL(int32 L, int32 R) noexcept {
      ASSERT(L > 0 && L <= X_ROWS);
      const ClusterRow &r = x_rows[L - 1];
      return (r.L == L) & (r.R == R);
    }
    
    [[nodiscard]] inline bool CLUSTONR(int32 L, int32 R) noexcept {
      ASSERT(R > 0 && R <= X_ROWS);
      const ClusterRow &r = x_rows[R - 1];
      return (r.L == L) & (r.R == R);
    }
    
    [[nodiscard]] inline bool ISCLUST(int32 L, int32 R) noexcept {
      // This function procedure returns value true if cluster <L,R> is in X;
      // otherwise it returns value false
      ASSERT(L > 0 && L <= X_ROWS);
      const ClusterRow &r_L = x_rows[L - 1];
      if ((r_L.L == L) & (r_L.R == R)) return true;
      
      ASSERT(L != R);
      ASSERT(R > 0 && R <= X_ROWS);
      const ClusterRow &r_R = x_rows[R - 1];
      return (r_R.L == L) & (r_R.R == R);
    }

    inline void CLEAR() noexcept {
      // Each cluster in X has an associated switch that is either cleared or
      // set.
      // This procedure clears every cluster switch in X.
      Xswitch.reset();
      xswitch_set_count = 0;
    }
        
    inline void SETSWX(std::size_t row) noexcept {
      // Only increment our counter on a 0 -> 1 transition
      ASSERT(row > 0 && row <= static_cast<size_t>(X_ROWS));
      if (!Xswitch[row]) {
        Xswitch.set(row, true);
        ++xswitch_set_count;
      }
    }

    [[nodiscard]] inline bool GETSWX(int32* row) noexcept {
      ASSERT(*row > 0 && *row <= X_ROWS);
      return Xswitch[*row];
    }

    [[nodiscard]] inline bool NOSWX(const std::size_t& n) noexcept {
      return xswitch_set_count == n;
    }
    
    inline void SETSW(int32 L, int32 R) noexcept {
      // If <L,R> is a cluster in X, 
      // this procedure sets the cluster switch for <L,R>.
      if (CLUSTONL(L, R)) {
        ++n_shared;
        SETSWX(L);
      } else if (CLUSTONR(L, R)) {
        ++n_shared;
        SETSWX(R);
      }
    }

    inline void UPDATE() noexcept {
      // This procedure inspects every cluster switch in X.
      // If the switch for cluster <L,R> is cleared, UPDATE deletes <L,R>
      // from X; thereafter ISCLUST(X,L,R) will return the value false.
      for (int32 i = X_ROWS; i--; ) {
        if (!(Xswitch[i])) {
          x_rows[i].L = 0;
          x_rows[i].R = 0;
        }
      }
    }

    [[nodiscard]] inline int32 SHARED() noexcept {
      // Used by COMCLUST in TreeDist::Day_1985.cpp
      return n_shared;
    }

    inline void ADDSHARED() noexcept {
      ++n_shared;
    }

    inline void XRESET() noexcept {
      // This procedure prepares X for an enumeration of its clusters
      enumeration = 0;
    }

    inline void NCLUS(int32* L, int32* R) noexcept {
      // This procedure returns the next cluster <L,R> in the current
      // enumeration of clusters in X.
      // If m clusters are in X, they are returned by the first m invocations
      // of NCLUS after initialization by XRESET; thereafter NCLUS returns the
      // invalid cluster <0,0>.
      ASSERT(enumeration < X_ROWS);
      const int32 row = static_cast<int32>(enumeration + 1); // rows are 1-based
      *L = X_left(row);
      *R = X_right(row);
      ++enumeration;
    }

  };

  inline ClusterTable::ClusterTable(Rcpp::List phylo) {

    const Rcpp::List rooted = TreeTools::root_on_node(phylo, 1);
    const Rcpp::IntegerMatrix edge = rooted["edge"];

    // BEGIN
    n_internal = rooted["Nnode"]; // = M
    Rcpp::CharacterVector leaf_labels = rooted["tip.label"];
    if (leaf_labels.length() > int(ct_max_leaves_heap)) {
      std::ostringstream msg;
      msg << "Tree has too many leaves (>" << ct_max_leaves_heap << "). "
          << "Contact the 'TreeTools' maintainer.";
      Rcpp::stop(msg.str());
    }
    ASSERT(ct_max_leaves_heap <= std::numeric_limits<int32>::max());
    n_leaves = int32(leaf_labels.length()); // = N
    if (double(edge.nrow()) > double(std::numeric_limits<int32>::max())) {
      std::ostringstream msg;
      msg << "Tree has too many edges (" << edge.nrow() << " > " << 
        std::numeric_limits<int32>::max() << "). " << 
          "Contact the 'TreeTools' maintainer.";
      Rcpp::stop(msg.str());
    }
    n_edge = int32(edge.nrow());
    const int32 n_vertex = M() + N();
    Tlen = 2 * n_vertex;
    Tlen_short = Tlen - (2 * 3);
    T = std::vector<int32>(Tlen);
    T_idx = 0;

    resize_uninitialized(leftmost_leaf, n_vertex);
    resize_uninitialized(visited_nth, n_leaves);
    resize_uninitialized(internal_label, 1 + n_leaves); // We're not using -1
    ASSERT(leftmost_leaf.size() >= static_cast<size_t>(n_vertex));
    ASSERT(visited_nth.size()   >= static_cast<size_t>(n_leaves));
    ASSERT(internal_label.size()>= static_cast<size_t>(1 + n_leaves));
    internal_label_ptr = internal_label.data();
    int32 n_visited = 0;
    std::vector<int32> weights(1 + n_vertex);

    for (int32 i = 1; i != n_leaves + 1; ++i) {
      SET_LEFTMOST(i, i); // Line 402
      weights[i] = 0;
    }
    for (int32 i = 1 + n_leaves; i != 1 + n_vertex; ++i) {
      SET_LEFTMOST(i, 0);
      weights[i] = 0;
    }
    for (int32 i = n_edge; i--; ) {
      const int32
        parent_i = int32(edge(i, 0)),
        child_i = int32(edge(i, 1))
      ;
      if (!GET_LEFTMOST(parent_i)) {
        SET_LEFTMOST(parent_i, GET_LEFTMOST(child_i));
      }
      if (is_leaf(child_i)) {
        VISIT_LEAF(&child_i, &n_visited);
        ++weights[parent_i];
        ENTER(child_i, 0);
      } else {
        weights[parent_i] += 1 + weights[child_i];
        ENTER(child_i, weights[child_i]);
      }
    }
    ENTER(int32(edge(0, 0)), weights[edge(0, 0)]);

    // BUILD Cluster table
    X_ROWS = n_leaves;
    x_rows = std::vector<ClusterRow>(X_ROWS);
    
    // Initialize Xswitch with appropriate size
    // Uses stack allocation for small trees, heap for large ones
    Xswitch.resize(n_leaves + 1);

    // This procedure constructs in X descriptions of the clusters in a
    // rooted tree described by the postorder sequence T with weights,
    // BUILD assigns each leaf an internal label so that every cluster
    // is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
    // simply described by a pair <L,R> of internal labels.

    TRESET();
    for (int32 i = 1; i != N(); ++i) {
      setX_left(i, 0);
      setX_right(i, 0);
    }
    int32 leafcode = 0, v, w, L, R = UNINIT, loc;

    NVERTEX(&v, &w);
    while (v) {
      if (is_leaf(v)) {
        ++leafcode;
        // We prepared the encoder in an earlier step, so need no X[v, 3] <- leafcode
        R = leafcode;
        NVERTEX(&v, &w);
      } else {
        L = ENCODE(LEFTLEAF());
        NVERTEX(&v, &w);
        loc = w == 0 ? R : L;
        setX_left(loc, L);
        setX_right(loc, R);
      }
    }
  }
}

#endif
