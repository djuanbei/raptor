#include "mcfcg.hpp"

namespace raptor {
namespace mcmcf {
KLUsolver::~KLUsolver() {
  if (!first) {
#ifdef USING_KLU
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);
#endif
  }
}

void KLUsolver::update(vector<SparseMatrixElem> &els, int n) {
#ifdef USING_KLU
  elements = els;
  if (!first) {
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);
  }
  first = false;
  Ap = new int[n + 1];
  Ai = new int[elements.size()];
  Ax = new double[elements.size()];
  stable_sort(elements.begin(), elements.end());

  int column = 0;
  Ap[column] = 0;
  int row = 0;
  double value = 0;
  int nz = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    int tempColumn = elements[i].column;
    int tempRow = elements[i].row;
    if (column == tempColumn && row == tempRow) {
      value += elements[i].value;
    } else {
      if (fabs(value) > 1e-6) {
        Ai[nz] = row;
        Ax[nz] = value;
        nz++;
      }

      if (column != tempColumn) {
        for (int j = column; j < tempColumn; j++) {
          Ap[j + 1] = nz;
        }
      }
      column = tempColumn;
      row = tempRow;
      value = elements[i].value;
    }
  }
  if (fabs(value) > 1e-6) {
    Ai[nz] = row;
    Ax[nz] = value;
    nz++;
  }

  for (int j = column; j < n; j++) {
    Ap[j + 1] = nz;
  }
  Symbolic = klu_analyze(n, Ap, Ai, &Common);
  Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);

  nonzeroNum = nz;
  dim = n;
  delete[] Ai;
  delete[] Ap;
  delete[] Ax;
#endif
}
bool KLUsolver::solve(double *b) {
#ifdef USING_KLU
  return klu_solve(Symbolic, Numeric, dim, 1, b, &Common) == 1;
#else
  return true;
#endif
}
bool KLUsolver::tsolve(double *b) {
#ifdef USING_KLU
  return klu_tsolve(Symbolic, Numeric, dim, 1, b, &Common) == 1;
#else
  return true;
#endif
}
}
}
