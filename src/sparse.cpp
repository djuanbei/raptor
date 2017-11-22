#include "sparse.h"
#include <cassert>
#include "klu.h"
#include "util.h"
namespace raptor {
namespace sparse {
using namespace std;
SparseSolver::SparseSolver(vector<sparseMatrixElem>& elements) {
  cs* T = cs_spalloc(0, 0, 1, 1, 1); /* allocate result */
  for (vector<sparseMatrixElem>::iterator it = elements.begin();
       it != elements.end(); it++) {
    cs_entry(T, (csi)it->row, (csi)it->column, it->value);
  }
  A = cs_compress(T);
  TA = cs_transpose(A, 1);
  cs_spfree(T); /* clear T */
}

SparseSolver::~SparseSolver() {
  cs_spfree(A);
  cs_spfree(TA);
}

void SparseSolver::minComputableProjection(const SparseVector& b, set<int>& I,
                                           set<int>& J) const {
  csi n = A->n, m = A->m, *Ap = A->p, *Ai = A->i, *tAp, *tAi;

  tAp = TA->p;
  tAi = TA->i;

  set<int> nI, nJ;
  I.insert(b.locs.begin(), b.locs.end());

  nI = I;
  while (!nI.empty()) {
    nJ.clear();
    for (set<int>::iterator it = nI.begin(); it != nI.end(); it++) {
      int start = tAp[*it];
      int end = tAp[*it + 1];
      for (int i = start; i < end; i++) {
        if (J.find(tAi[i]) == J.end()) {
          nJ.insert(tAi[i]);
        }
      }
    }

    if (nJ.empty()) {
      break;
    }
    J.insert(nJ.begin(), nJ.end());
    nI.clear();
    for (set<int>::iterator it = nJ.begin(); it != nJ.end(); it++) {
      int start = Ap[*it];
      int end = Ap[*it + 1];
      for (int i = start; i < end; i++) {
        if (I.find(Ai[i]) == I.end()) {
          nI.insert(Ai[i]);
        }
      }
    }
    I.insert(nI.begin(), nI.end());
  }
}

bool SparseSolver::locSolver(SparseVector& sB) const {
  csi *Ap = A->p, *Ai = A->i;
  double* X = A->x;
  set<int> I, J;

  minComputableProjection(sB, I, J);

  assert(I.size() >= J.size());

  vector<int> vecI(I.begin(), I.end());
  vector<int> vecJ(J.begin(), J.end());
  int smallest = vecI.front();
  int bigest = vecI.back();

  cs* B = cs_spalloc(0, 0, 1, 1, 1);

  for (size_t j = 0; j < vecJ.size(); j++) {
    int k = vecJ[j];
    int start = Ap[k];
    int end = Ap[k + 1];
    for (int l = start; l < end; l++) {
      int p = Ai[l];
      if (p >= smallest && p <= bigest) {
        int i = binfind(vecI.begin(), vecI.end(), p);
        if (i > -1) {
          cs_entry(B, i, j, X[l]);
        }
      }
    }
  }

  double* y = new double[J.size()];
  fill(y, y + J.size(), 0);
  for(size_t j=0; j< sB.locs.size(); j++){
    int k = sB.locs[j];
    int i = binfind(vecI.begin(), vecI.end(), k);
    assert(i > -1);
    y[i] =sB.values[j];
  }


  cs* M = cs_compress(B);

  cs_spfree(B);

  if (vecI.size() != vecJ.size()) {
    cs* TB = cs_transpose(M, 1);
    cs* B = cs_multiply(M, TB);
    cs_spfree(TB);
    cs_spfree(M);
    M = B;
  }

  klu_symbolic* Symbolic;
  klu_numeric* Numeric;
  klu_common Common;
  klu_defaults(&Common);
  int n = M->n;
  int* p = new int[n + 1];
  for (int i = 0; i <= n; i++) {
    p[i] = M->p[i];
  }
  int num = M->p[n];
  int* index = new int[num];
  for (int i = 0; i < num; i++) {
    index[i] = M->i[i];
  }

  Symbolic = klu_analyze(M->n, p, index, &Common);
  Numeric = klu_factor(p, index, M->x, Symbolic, &Common);

  int re = klu_solve(Symbolic, Numeric, M->n, 1, y, &Common);

  klu_free_symbolic(&Symbolic, &Common);
  klu_free_numeric(&Numeric, &Common);

  cs_spfree(M);
  if (1 != re) {
    delete[] y;
    delete[] p;
    delete[] index;
    return false;
  }
  /**
   * lifting
   *
   */

  sB.clear();
  for (size_t i = 0; i < vecJ.size(); i++) {
    if(fabs(y[i])>1e-6){
      int k = vecJ[i];
      sB.locs.push_back(k);
      sB.values.push_back(y[i]);
    }
  }
  delete[] y;
  delete[] p;
  delete[] index;

  return true;
}
bool SparseSolver::locSolver(double* b) const {

  SparseVector sB;
  for (int i = 0; i < A->n; i++) {
    if (fabs(b[i]) > 1e-6) {
      sB.locs.push_back(i);
      sB.values.push_back(b[i]);
    }
  }
  if(! locSolver(sB)){
    return false;
  }
  fill(b, b + A->n, 0);
  for(size_t i=0; i< sB.locs.size(); i++){
    b[sB.locs[i]]=sB.values[i];
  }
  return true;
}

}
}
