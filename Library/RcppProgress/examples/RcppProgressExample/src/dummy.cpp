/* regression test
 dummy source to test impact of multiple inclusions
 of RcppProgress header files.
*/


#include "progress.hpp"

void toto() {
  Progress p(100, false);
}

