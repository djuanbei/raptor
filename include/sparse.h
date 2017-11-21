/**
 * @file   sparse.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Mon Nov 20 19:36:04 2017
 *
 * @brief  sparse linear equation  solver Ax=b where A is a sparse matrix and b
 * is sparse vector
 *
 *
 */

#include "config.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "csparse/cs.h"
#ifdef __cplusplus
}
#endif
#include <set>
#include <vector>

namespace raptor {
namespace sparse {
using namespace std;

class SparseSolver {
 private:
  cs *A;
  void minComputableProjection(const vector<pair<int, double> >& b, set<int>& I,
                               set<int>& J) const;

 public:
  SparseSolver(vector<sparseMatrixElem>& elements);
  ~SparseSolver();
  bool locSolver(double *b) const;
  bool incSolver(double *b) const;
};
}
}
