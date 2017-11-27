#include "sparse.h"
#include <cassert>
#include "graph.h"
#include "graphalg.hpp"
#include "klu.h"
#include "util.h"
namespace raptor {
namespace sparse {
using namespace std;

subMatrix subMatrix::operator+(const subMatrix& other) const {
  subMatrix re = other;
  re.add(other);

  return re;
}

void subMatrix::add(const subMatrix& other) {
  rows.insert(rows.end(), rows.begin(), rows.end());
  columns.insert(columns.end(), columns.begin(), columns.end());
  sort(rows.begin(), rows.end());
  sort(columns.begin(), columns.end());

#ifdef DEBUG
  for (size_t i = 0; i + 1 < rows.size(); i++) {
    assert(rows[i] != rows[i + 1]);
  }
  for (size_t i = 0; i + 1 < columns.size(); i++) {
    assert(columns[i] != columns[i + 1]);
  }
#endif
}

SparseSolver::SparseSolver() { A = TA = NULL; }

void SparseSolver::update(const vector<sparseMatrixElem>& elements) {
  
  cs* T = cs_spalloc(0, 0, 1, 1, 1); /* allocate result */
  for (vector<sparseMatrixElem>::const_iterator it = elements.begin();
       it != elements.end(); it++) {
    cs_entry(T, (csi)it->row, (csi)it->column, it->value);
  }
  if(NULL!=A){
    cs_spfree(A);
  }
  if(NULL!=TA){
    cs_spfree(TA);
  }
  
  A = cs_compress(T);
  TA = cs_transpose(A, 1);
  cs_spfree(T); /* clear T */
  dim = A->n;

  computableConnectedComponent();
}

SparseSolver::~SparseSolver() {
  if (NULL != A) {
    cs_spfree(A);
  }
  if (NULL != TA) {
    cs_spfree(TA);
  }
}

void SparseSolver::minComputableProjection(const SparseVector& b,
                                           vector<int>& I,
                                           vector<int>& J) const {
  I.clear();
  J.clear();
  set<int> passSubs;
  for (vector<int>::const_iterator it = b.locs.begin(); it != b.locs.end();
       it++) {
    passSubs.insert(subIndex[*it]);
  }

  for (set<int>::iterator it = passSubs.begin(); it != passSubs.end(); it++) {
    I.insert(I.end(), subs[*it].rows.begin(), subs[*it].rows.end());
    J.insert(J.end(), subs[*it].columns.begin(), subs[*it].columns.end());
  }

  for (vector<int>::iterator it = J.begin(); it != J.end(); it++) {
    *it -= dim;
  }
  sort(I.begin(), I.end());
  sort(J.begin(), J.end());
}

void SparseSolver::tminComputableProjection(const SparseVector& b,
                                            vector<int>& I,
                                            vector<int>& J) const {
  I.clear();
  J.clear();
  set<int> passSubs;
  for (vector<int>::const_iterator it = b.locs.begin(); it != b.locs.end();
       it++) {
    passSubs.insert(subIndex[*it + dim]);
  }

  for (set<int>::iterator it = passSubs.begin(); it != passSubs.end(); it++) {
    I.insert(I.end(), subs[*it].columns.begin(), subs[*it].columns.end());

    J.insert(J.end(), subs[*it].rows.begin(), subs[*it].rows.end());
  }

  for (vector<int>::iterator it = J.begin(); it != J.end(); it++) {
    *it -= dim;
  }
  sort(I.begin(), I.end());
  sort(J.begin(), J.end());
}
bool SparseSolver::locSolver(SparseVector& sB, csi* Ap, csi* Ai, double* X,
                             vector<int>& vecI, vector<int>& vecJ) const {
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

  double* y = new double[vecJ.size()];
  fill(y, y + vecJ.size(), 0);
  for (size_t j = 0; j < sB.locs.size(); j++) {
    int k = sB.locs[j];
    int i = binfind(vecI.begin(), vecI.end(), k);
    assert(i > -1);
    y[i] = sB.values[j];
  }

  cs* M = cs_compress(B);

  cs_spfree(B);

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
    if (fabs(y[i]) > 1e-6) {
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

void SparseSolver::computableConnectedComponent() {
  if (NULL != TA) {
    subs.clear();

    csi n = dim, *Ap = TA->p, *Ai = TA->i;
    subIndex.resize(2 * n);
    fill(subIndex.begin(), subIndex.end(), 0);
    vector<int> srcs, snks;

    for (csi i = 0; i < n; i++) {
      int start = Ap[i];
      int end = Ap[i + 1];
      for (int j = start; j < end; j++) {
        int k = Ai[j] + n;
        srcs.push_back(i);
        snks.push_back(k);

        srcs.push_back(k);
        snks.push_back(i);
      }
    }

    simple_graph graph;
    graph.initial(srcs, snks);
    vector<char> pass(n, 0);
    for (int i = 0; i < n; i++) {
      if (!pass[i]) {
        subMatrix temp;
        vector<int> passNodes;
        bfs_search(graph, i, passNodes);
        sort(passNodes.begin(), passNodes.end());
        for (vector<int>::iterator it = passNodes.begin();
             it != passNodes.end(); it++) {
          subIndex[*it] = subs.size();
        }

        vector<int>::iterator it =
            lower_bound(passNodes.begin(), passNodes.end(), n);
        temp.rows.insert(temp.rows.end(), passNodes.begin(), it);
        for (vector<int>::iterator it = temp.rows.begin();
             it != temp.rows.end(); it++) {
          pass[*it] = 1;
        }
        temp.columns.insert(temp.columns.end(), it, passNodes.end());

        subs.push_back(temp);
      }
    }
  }
}

bool SparseSolver::locSolver(SparseVector& sB) const {
  csi *Ap = A->p, *Ai = A->i;
  double* X = A->x;
  vector<int> vecI, vecJ;

  minComputableProjection(sB, vecI, vecJ);

  assert(vecI.size() == vecJ.size());

  return locSolver(sB, Ap, Ai, X, vecI, vecJ);
}
bool SparseSolver::locSolver(double* b) const {
  SparseVector sB;
  for (int i = 0; i < dim; i++) {
    if (fabs(b[i]) > 1e-6) {
      sB.locs.push_back(i);
      sB.values.push_back(b[i]);
    }
  }
  if (!locSolver(sB)) {
    return false;
  }
  fill(b, b + dim, 0);
  for (size_t i = 0; i < sB.locs.size(); i++) {
    b[sB.locs[i]] = sB.values[i];
  }
  return true;
}

bool SparseSolver::tlocSolver(SparseVector& sB) const {
  csi *Ap = TA->p, *Ai = TA->i;
  double* X = TA->x;
  vector<int> vecI, vecJ;

  tminComputableProjection(sB, vecI, vecJ);

  assert(vecI.size() == vecJ.size());
  return locSolver(sB, Ap, Ai, X, vecI, vecJ);
}

bool SparseSolver::tlocSolver(double* b) const {
  SparseVector sB;
  for (int i = 0; i < dim; i++) {
    if (fabs(b[i]) > 1e-6) {
      sB.locs.push_back(i);
      sB.values.push_back(b[i]);
    }
  }
  if (!tlocSolver(sB)) {
    return false;
  }
  fill(b, b + dim, 0);
  for (size_t i = 0; i < sB.locs.size(); i++) {
    b[sB.locs[i]] = sB.values[i];
  }
  return true;
}

bool SparseSolver::incSolver(const double *initSolution, double *b) const{
  csi n = dim, *Ap = TA->p, *Ai = TA->i;

  double *X=TA->x;
  for (csi i = 0; i < n; i++) {

    int start = Ap[i];
    int end = Ap[i + 1];
    for (int k = start; k < end; k++) {
      int j = Ai[k];
      b[i]-=X[k]*initSolution[j];
    }
  }

  if(!locSolver(b)){
    return false;
  }
  for(int i=0; i< dim; i++){
    b[i]+=initSolution[i];
  }
  return true;
  
}


bool SparseSolver::tincSolver(const double *initSolution, double *b) const{
  csi n = dim, *Ap = A->p, *Ai = A->i;

  double *X=A->x;
  for (csi i = 0; i < n; i++) {

    int start = Ap[i];
    int end = Ap[i + 1];
    for (int k = start; k < end; k++) {
      int j = Ai[k];
      b[i]-=X[k]*initSolution[j];
    }
  }

  if(!tlocSolver(b)){
    return false;
  }
  for(int i=0; i< dim; i++){
    b[i]+=initSolution[i];
  }
  return true;
  
}


}
}
