#ifndef TreeTools_assert_
#define TreeTools_assert_
#include <sstream>

#if defined(DEBUG) && !defined(PROFILE)
#define ASSERT(x) if (!(x)) {                                  \
std::ostringstream oss;                                        \
oss << "Assertion failed: " << #x;                             \
Rcpp::stop(oss.str());                                         \
}
#else
#define ASSERT(x) ((void)0)
#endif

#endif
