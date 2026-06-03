#ifndef TreeTools_tree_number_
#define TreeTools_tree_number_

// Unique integer identifiers for unrooted bifurcating tree topologies.
// Uses a 256-bit unsigned integer to support up to 51 leaves.
//
// Based on the mixed-base encoding of John Tromp (Li et al. 1996):
// each n-leaf tree is uniquely identified by a non-negative integer
// computed from the sequence of edge insertions that build the tree.
//
// The encoding/decoding algorithms are adapted from TreeTools' int_to_tree.cpp.

#include <cstdint>
#include <functional>
#include "types.h"

namespace TreeTools {

constexpr intx TREE_NUM_MAX_TIP = 51;
constexpr intx TREE_NUM_MAX_NODE = TREE_NUM_MAX_TIP * 2 - 1;

// ============================================================
// tree_num_t: lightweight 256-bit unsigned integer
// ============================================================
// 4 × uint64_t, little-endian (w[0] = least significant word).
// Only implements operations needed for tree number encoding/decoding.

struct tree_num_t {
  uint64_t w[4];

  tree_num_t() : w{0, 0, 0, 0} {}
  explicit tree_num_t(uint64_t v) : w{v, 0, 0, 0} {}

  bool operator==(const tree_num_t& o) const {
    return w[0] == o.w[0] && w[1] == o.w[1] &&
           w[2] == o.w[2] && w[3] == o.w[3];
  }
  bool operator!=(const tree_num_t& o) const { return !(*this == o); }

  // Full 256-bit addition.
  tree_num_t& operator+=(const tree_num_t& o) {
    uint64_t carry = 0;
    for (int i = 0; i < 4; ++i) {
      uint64_t a = w[i], b = o.w[i];
      uint64_t sum = a + b;
      uint64_t c1 = (sum < a) ? 1ULL : 0ULL;
      uint64_t sum2 = sum + carry;
      uint64_t c2 = (sum2 < sum) ? 1ULL : 0ULL;
      w[i] = sum2;
      carry = c1 + c2;
    }
    return *this;
  }

  // Multiply *this by v (v < 2^64). Returns *this.
  tree_num_t& mul_small(uint64_t v) {
#ifdef __SIZEOF_INT128__
    __uint128_t carry = 0;
    for (int i = 0; i < 4; ++i) {
      __uint128_t prod = static_cast<__uint128_t>(w[i]) * v + carry;
      w[i] = static_cast<uint64_t>(prod);
      carry = prod >> 64;
    }
#else
    // Fallback without 128-bit support: process 32-bit halves.
    // Correct when v < 2^32 (always true for tree bases, max = 2*51-3 = 99).
    uint64_t carry = 0;
    for (int i = 0; i < 4; ++i) {
      uint64_t lo32 = w[i] & 0xFFFFFFFFULL;
      uint64_t hi32 = w[i] >> 32;
      uint64_t prod_lo = lo32 * v + carry;
      uint64_t prod_hi = hi32 * v + (prod_lo >> 32);
      w[i] = (prod_hi << 32) | (prod_lo & 0xFFFFFFFFULL);
      carry = prod_hi >> 32;
    }
#endif
    return *this;
  }

  // Add v (< 2^64) to *this. Returns *this.
  tree_num_t& add_small(uint64_t v) {
    uint64_t sum = w[0] + v;
    uint64_t carry = (sum < w[0]) ? 1ULL : 0ULL;
    w[0] = sum;
    for (int i = 1; i < 4 && carry; ++i) {
      sum = w[i] + carry;
      carry = (sum < w[i]) ? 1ULL : 0ULL;
      w[i] = sum;
    }
    return *this;
  }

  // Divide *this by v (v > 0), return remainder. Modifies *this to quotient.
  uint64_t divmod_small(uint64_t v) {
#ifdef __SIZEOF_INT128__
    __uint128_t rem = 0;
    for (int i = 3; i >= 0; --i) {
      __uint128_t cur = (rem << 64) | w[i];
      w[i] = static_cast<uint64_t>(cur / v);
      rem = cur % v;
    }
    return static_cast<uint64_t>(rem);
#else
    // Fallback: chain 32-bit divisions (correct for v < 2^32).
    uint64_t rem = 0;
    for (int i = 3; i >= 0; --i) {
      uint64_t hi = (rem << 32) | (w[i] >> 32);
      uint64_t q_hi = hi / v;
      rem = hi % v;
      uint64_t lo = (rem << 32) | (w[i] & 0xFFFFFFFFULL);
      uint64_t q_lo = lo / v;
      rem = lo % v;
      w[i] = (q_hi << 32) | q_lo;
    }
    return rem;
#endif
  }
};

// Hash functor for unordered containers.
struct tree_num_hash {
  std::size_t operator()(const tree_num_t& v) const {
    std::size_t h = std::hash<uint64_t>{}(v.w[0]);
    h ^= std::hash<uint64_t>{}(v.w[1]) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<uint64_t>{}(v.w[2]) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<uint64_t>{}(v.w[3]) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

// ============================================================
// Encoding: postorder edge arrays → tree_num_t
// ============================================================
// parent[] and child[] are 1-indexed, in postorder with consecutive
// edge pairs sharing the same parent. Tree rooted on tip 1.
// n_tip must be >= 2 and <= TREE_NUM_MAX_TIP.

inline tree_num_t edges_to_tree_number(
    const intx* parent, const intx* child, intx n_tip) {
  const intx all_node = n_tip + n_tip - 1;
  const intx n_edge = n_tip + n_tip - 2;

  intx smallest_below[TREE_NUM_MAX_NODE];
  intx parent_of[TREE_NUM_MAX_NODE];
  intx prime_id[TREE_NUM_MAX_NODE];
  intx index[TREE_NUM_MAX_TIP];

  for (intx i = 0; i < all_node; ++i) {
    smallest_below[i] = i;
    prime_id[i] = i;
  }

  for (intx i = 0; i < n_edge - 2; i += 2) {
    const intx nd = parent[i] - 1;
    const intx lc = child[i] - 1;
    const intx rc = child[i + 1] - 1;

    smallest_below[nd] = (smallest_below[rc] < smallest_below[lc])
      ? smallest_below[rc] : smallest_below[lc];
    prime_id[nd] = (smallest_below[lc] > smallest_below[rc])
      ? smallest_below[lc] : smallest_below[rc];
    parent_of[lc] = parent_of[rc] = nd;

    for (intx at = smallest_below[nd]; at != nd; at = parent_of[at]) {
      if (prime_id[at] < prime_id[nd]) {
        index[prime_id[nd]] = prime_id[at] + (at < n_tip ? 0 : n_tip);
      }
    }
  }

  tree_num_t num;
  tree_num_t multiplier(1);
  for (intx i = 3; i != n_tip; ++i) {
    intx ie = index[i];
    if (ie < n_tip) {
      --ie;
    } else {
      ie += i - (n_tip + 3);
    }
    tree_num_t term(multiplier);
    term.mul_small(static_cast<uint64_t>(ie));
    num += term;
    multiplier.mul_small(static_cast<uint64_t>(i + i - 3));
  }
  return num;
}

// ============================================================
// Decoding: tree_num_t → parent vector
// ============================================================
// Writes parent_out[0..2*n_tip-3], where parent_out[i] is the
// 1-indexed parent of node (i+1). Root node = 2*n_tip - 1.

inline void tree_number_to_parent(
    tree_num_t tree_id, intx n_tip, intx* parent_out) {
  const intx root_node = n_tip + n_tip - 1;
  const intx prime = n_tip - 2;

  parent_out[0] = root_node;
  parent_out[1] = root_node;

  for (intx i = 2; i != n_tip; ++i) {
    const intx base = i + i - 3;
    const intx i_prime = i + prime;

    intx where = static_cast<intx>(
      tree_id.divmod_small(static_cast<uint64_t>(base))) + 1;
    if (where >= i) {
      where += prime + 2 - i;
    }

    parent_out[i_prime] = parent_out[where];
    parent_out[i] = i_prime + 1;
    parent_out[where] = i_prime + 1;
  }
}

} // namespace TreeTools

#endif
