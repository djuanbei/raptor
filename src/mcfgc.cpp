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
  delete[] Ai;
  delete[] Ap;
  delete[] Ax;
}

void KLUsolver::update(vector<SparseMatrixElem> &elements, int n) {
#ifdef USING_KLU
  if(n>cap|| (int)elements.size()>2*cap){
    int d=max(n,(int)elements.size()/2);
    reScale(d);
  }
  if (!first) {
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);
  }
  first = false;

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
void KLUsolver::reScale(int s){
  cap=s;
  cap+=2;
  cap*=1.3;
  if(NULL!=Ap){
    delete[] Ap;
    delete[] Ai;
    delete[] Ax;
  }
  Ap=new int[cap+1];
  Ai=new int[2*cap];
  Ax=new double[2*cap];
}
}
}
