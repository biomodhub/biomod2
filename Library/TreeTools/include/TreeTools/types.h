#ifndef TreeTools_types_
#define TreeTools_types_

#include <cstdint> /* for uint_fast32_t, int_fast32_t, INTPTR_MAX */
#include <limits> /* for numeric_limits */

using int16 = int_fast16_t;
using int32 = int_fast32_t;
using intx = int_fast32_t;
using uintx = uint_fast32_t;

namespace TreeTools {
  const intx INTX_MAX = intx(std::numeric_limits<intx>::max());
  const uintx UINTX_MAX = uintx(std::numeric_limits<uintx>::max());
}

#endif
