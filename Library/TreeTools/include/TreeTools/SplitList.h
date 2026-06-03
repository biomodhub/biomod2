#ifndef _TREETOOLS_SPLITLIST_H
#define _TREETOOLS_SPLITLIST_H

#include <Rcpp/Lightest>
#include <stdexcept> /* for errors */
#include <vector>    /* for heap allocation */
#include <algorithm> /* for std::fill */

#include "assert.h" /* for ASSERT */
#include "types.h" /* for int16, int32 */

using splitbit = uint_fast64_t;

#define R_BIN_SIZE int16(8)
#define SL_BIN_SIZE int16(64)
#define SL_MAX_BINS int16(511)

// Keep SL_MAX_TIPS + 2 within int16 range (32767) for downstream packages
#define SL_MAX_TIPS (SL_BIN_SIZE * SL_MAX_BINS) // 32704
#define SL_MAX_SPLITS (SL_MAX_TIPS - 3) 

/* Stack allocation thresholds.
 * Trees with n_splits <= SL_STACK_SPLITS AND n_bins <= SL_STACK_BINS
 * use fast stack arrays; larger trees fall back to heap allocation.
 * Kept at the pre-v1.16 values to avoid bloating SplitList objects.
 */
#define SL_STACK_BINS int16(32)
#define SL_STACK_SPLITS int16(SL_BIN_SIZE * SL_STACK_BINS - 3) // 2045

#define INLASTBIN(n, size) int16((size) - int16((size) - int16((n) % (size))) % (size))
#define INSUBBIN(bin, offset)                                  \
  splitbit(x(split, ((bin) * input_bins_per_bin) + (offset)))
#define INBIN(r_bin, bin) ((INSUBBIN((bin), (r_bin))) << (R_BIN_SIZE * (r_bin)))

namespace TreeTools {

  constexpr int input_bins_per_bin = SL_BIN_SIZE / R_BIN_SIZE;

  template<typename T>
  [[nodiscard]] constexpr T power_of_two(int bit_pos) noexcept {
    static_assert(std::is_unsigned_v<T>, "Use unsigned types for bit operations");
    assert(bit_pos >= 0 && bit_pos < int(sizeof(T) * 8));
    return T(1) << bit_pos;
  }
  
// Hardware POPCNT: available on all x86-64 since 2008 (Nehalem / Barcelona).
  // Inline asm emits the instruction directly, without requiring -mpopcnt.
#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__)
  inline int32 count_bits(splitbit x) {
    uint64_t result;
    __asm__ ("popcnt %1, %0" : "=r" (result) : "r" (x));
    return static_cast<int32>(result);
  }
#elif defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
  inline int32 count_bits(splitbit x) {
    return static_cast<int32>(__popcnt64(x));
  }
#elif defined(__GNUC__) || defined(__clang__)
  // Non-x86 (ARM, etc.): builtin maps to efficient native instruction
  inline int32 count_bits(splitbit x) {
    return static_cast<int32>(__builtin_popcountll(x));
  }
#else
  inline int32_t count_bits(splitbit x) {
    int32_t count = 0;
    while (x != 0) {
      x &= (x - 1);
      ++count;
    }
    return count;
  }
#endif

  class SplitList {
  public:
    int32 n_splits;
    int32 n_bins;
    int32* in_split;
    splitbit** state;
  
  private:
    /* STACK STORAGE (Fast path for small trees ≤ SL_STACK_SPLITS splits) */
    int32 stack_in_split[SL_STACK_SPLITS];
    splitbit stack_state[SL_STACK_SPLITS][SL_STACK_BINS];
    splitbit* stack_rows[SL_STACK_SPLITS];
    
    /* HEAP STORAGE (Large trees) */
    std::vector<int32> heap_in_split;
    std::vector<splitbit> heap_data;      
    std::vector<splitbit*> heap_rows;     
    
  public:
    SplitList(const Rcpp::RawMatrix &x) {
      
      const double n_rows = static_cast<double>(x.rows());
      
      /* Check limits */
      if (n_rows > static_cast<double>(std::numeric_limits<int32>::max())) {
        Rcpp::stop("Too many splits (exceeds int32 limit).");                   // #nocov
      }
      
      n_splits = int32(x.rows());
      ASSERT(n_splits >= 0);
      
      const int32 n_input_bins = int32(x.cols());
      ASSERT(n_input_bins > 0);
      n_bins = int32(n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;
      
      bool use_heap = (n_splits > SL_STACK_SPLITS) || (n_bins > SL_STACK_BINS);
      
      if (use_heap) {
        heap_in_split.resize(n_splits, 0);
        in_split = heap_in_split.data();
        
        size_t total_elements = static_cast<size_t>(n_splits) *
          static_cast<size_t>(n_bins);
        heap_data.resize(total_elements); 
        
        heap_rows.resize(n_splits);
        for (int32 i = 0; i < n_splits; ++i) {
          heap_rows[i] = &heap_data[i * n_bins];
        }
        state = heap_rows.data();
        
      } else {
        in_split = stack_in_split;
        
        for (int32 i = 0; i < n_splits; ++i) {
          stack_rows[i] = stack_state[i];
          in_split[i] = 0;
        }
        state = stack_rows;
      }
      
      
      for (int32 bin = 0; bin < n_bins - 1; ++bin) {
        const int32 bin_offset = bin * input_bins_per_bin;
        
        for (int32 split = 0; split < n_splits; ++split) {
          splitbit combined = splitbit(x(split, bin_offset));
          
          for (int32 input_bin = 1; input_bin < input_bins_per_bin; ++input_bin) {
            combined |= splitbit(x(split, bin_offset + input_bin)) <<
              (R_BIN_SIZE * input_bin);
          }
          
          state[split][bin] = combined;
          in_split[split] += count_bits(combined);
        }
      }
      
      const int32 last_bin = n_bins - 1;
      const int32 raggedy_bins = INLASTBIN(n_input_bins, R_BIN_SIZE);
      
      for (int32 split = 0; split < n_splits; ++split) {
        state[split][last_bin] = INSUBBIN(last_bin, 0);
        
        for (int32 input_bin = 1; input_bin < raggedy_bins; ++input_bin) {
          state[split][last_bin] += INBIN(input_bin, last_bin);
        }
        
        in_split[split] += count_bits(state[split][last_bin]);
      }
    }
    
    // Default destructor handles vector cleanup automatically
    ~SplitList() = default;
    
    // Disable copy/move to prevent pointer invalidation issues
    SplitList(const SplitList&) = delete;
    SplitList& operator=(const SplitList&) = delete;
  };
}

#endif
