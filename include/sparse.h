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

struct SparseVector {
  vector<int> locs;
  vector<double> values;

  void clear(){
    locs.clear();
    values.clear();
  }
};

class SparseSolver {
 private:
  cs *lastA;
  cs *lastTA;
  vector<int> lastRowIndex;
  vector<double> lastRightHandSide;

  cs *A;
  cs *TA;
  vector<int> rowIndex;
  vector<double> rightHandSide;

  void minComputableProjection(const SparseVector &b, set<int> &I,
                               set<int> &J) const;

 public:
  SparseSolver(vector<sparseMatrixElem> &elements);
  ~SparseSolver();

  /**
   * @brief local linear equation system solver
   *
   * @param b right hand side
   *
   * @return  rewrite solution to b
   */
  bool locSolver(SparseVector &b) const;

  /**
   * @brief local linear equation system solver
   *
   * @param b right hand side
   *
   * @return  rewrite solution to b
   */
  bool locSolver(double *b) const;

  /**
   * @brief local linear equation system solver for A^Tx=b
   *
   * @param b right hand side
   *
   * @return  rewrite solution to b
   */
  bool tlocSolver(SparseVector &b) const;

  /**
   * @brief local linear equation system solver for A^Tx=b
   *
   * @param b right hand side
   *
   * @return  rewrite solution to b
   */
  bool tlocSolver(double *b) const;

  /**
   * @brief incremental linear equation solver
   *
   * @param b  right hand side
   *
   * @return rewrite solution to b
   */
  bool incSolver(double *b) const;
};
}
}
