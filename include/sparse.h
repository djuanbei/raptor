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

#ifndef __SPARSE_LES_H
#define __SPARSE_LES_H

#include "config.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "csparse/cs.h"
#ifdef __cplusplus
}
#endif
#include <algorithm>
#include <cassert>
#include <set>
#include <vector>

namespace raptor {
namespace sparse {
using namespace std;

struct SparseVector {
  vector<int> locs;
  vector<double> values;

  void clear() {
    locs.clear();
    values.clear();
  }
};

struct subMatrix {
  vector<int> rows, columns;

  subMatrix operator+(const subMatrix &other) const;
  void add(const subMatrix &other);
  void update() {
    sort(rows.begin(), rows.end());
    sort(columns.begin(), columns.end());
  }
};

class SparseSolver {
 private:
  int dim;

  cs *A;
  cs *TA;

  vector<subMatrix> subs;
  vector<int> subIndex;

  void minComputableProjection(const SparseVector &b, vector<int> &I,
                               vector<int> &J) const;

  void tminComputableProjection(const SparseVector &b, vector<int> &I,
                                vector<int> &J) const;

  bool locSolver(SparseVector &sB, csi *Ap, csi *Ai, double *X,
                 vector<int> &vecI, vector<int> &vecJ) const;

  void computableConnectedComponent();

 public:
  SparseSolver();

  // SparseSolver(vector<SparseMatrixElem> &elements);

  ~SparseSolver();

  /**
   *
   * @brief update next matrix
   * @param elements
   *
   * @return
   */
  void update( vector<SparseMatrixElem> &elements);

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
  bool incSolver(const double *initSolution, double *b) const;

  /**
   * @brief incremental linear equation solver A^Tx=b
   *
   * @param b  right hand side
   *
   * @return rewrite solution to b
   */
  bool tincSolver(const double *initSolution, double *b) const;
};
}
}
#endif
