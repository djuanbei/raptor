
/**
 * @file   config.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Thu Nov  9 10:28:03 2017
 *
 * @brief  some util data structure
 *
 *
 */

#ifndef __CONFIG_H
#define __CONFIG_H
#include <limits>

#ifdef USING_LAPACK
/* DGESV prototype */
extern "C" {
void dgesv_(int *N, int *nrhs, double *a, int *lda, int *ipiv, double *b,
            int *ldb, int *info);

void sgesv_(int *N, int *nrhs, float *a, int *lda, int *ipiv, float *b,
            int *ldb, int *info);

// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
//// generate inverse of a matrix given its LU decomposition
// void dgetri_( int* N, double* A, int* lda, int* IPIV, double* WORK, int*
// lwork, int* INFO);
void dgetrs_(char *C, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
             double *B, int *LDB, int *INFO);

// LU decomoposition of a general matrix
void sgetrf_(int *M, int *N, float *A, int *lda, int *IPIV, int *INFO);
//// generate inverse of a matrix given its LU decomposition
// void dgetri_( int* N, double* A, int* lda, int* IPIV, double* WORK, int*
// lwork, int* INFO);
void sgetrs_(char *C, int *N, int *NRHS, float *A, int *LDA, int *IPIV,
             float *B, int *LDB, int *INFO);
}
#else

extern "C" {
static void dgesv_(int *N, int *nrhs, double *a, int *lda, int *ipiv, double *b,
                   int *ldb, int *info) {}

static void sgesv_(int *N, int *nrhs, float *a, int *lda, int *ipiv, float *b,
                   int *ldb, int *info) {}

// LU decomoposition of a general matrix
static void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO) {
}
//// generate inverse of a matrix given its LU decomposition
// void dgetri_( int* N, double* A, int* lda, int* IPIV, double* WORK, int*
// lwork, int* INFO);
static void dgetrs_(char *C, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
                    double *B, int *LDB, int *INFO) {}

// LU decomoposition of a general matrix
static void sgetrf_(int *M, int *N, float *A, int *lda, int *IPIV, int *INFO) {}
//// generate inverse of a matrix given its LU decomposition
// void dgetri_( int* N, double* A, int* lda, int* IPIV, double* WORK, int*
// lwork, int* INFO);
static void sgetrs_(char *C, int *N, int *NRHS, float *A, int *LDA, int *IPIV,
                    float *B, int *LDB, int *INFO) {}
}

#endif

#define callTime(t, exper)                                \
  auto t0 = std::chrono::system_clock::now();             \
  exper;                                                  \
  auto t1 = std::chrono::system_clock::now();             \
  std::chrono::duration<double, std::milli> ms = t1 - t0; \
  t = ms.count();

namespace raptor {
struct SparseMatrixElem {
  int column, row;
  double value;
  SparseMatrixElem() : column(0), row(0), value(0) {}
  bool operator<(const SparseMatrixElem &other) const {
    if (column < other.column) {
      return true;
    }
    if (column > other.column) {
      return false;
    }
    return row <= other.row;
  }
};

template <typename W>
W getInf(W) {
  return std::numeric_limits<W>::max() / 30.0;
}

enum LU_SOLVER {
  KLU = 0,
  LAPACK = 1,
  SPARSE=2
};

template <typename C>
struct Demand {
  int src, snk;
  C bandwidth;
  Demand() : src(0), snk(0), bandwidth(0) {}
};

enum ENTER_BASE_TYPE {
  PATH_T = 0,
  LINK_T = 1

};

enum EXIT_BASE_TYPE {
  DEMAND_T = 0,
  STATUS_LINK = 1,
  OTHER_LINK = 2

};

enum PivotSaturateMatrix {
  NOCHANGE = 0,
  ADDLINK = 1,
  DELETELINK = 2

};

struct Statistics_data {
  int iterator_num;
  int empty_iterator_num;
  double bestUpperobj;
  double estimee_opt_diff;
  int nzn;
  int snzn;
  double minLU;
  ENTER_BASE_TYPE etype;
  int enter;
  EXIT_BASE_TYPE exitt;
  int leaving;
  double start_time;
  double using_system_time;
  double objSpeed;
  double lpsolvertime;
  double shortestpathtime;
  double totalStaturateLink;
  double totalNonzero;
  PivotSaturateMatrix pivotType;
  int linkId;
  Statistics_data()
      : iterator_num(0),
        empty_iterator_num(0),
        estimee_opt_diff(std::numeric_limits<double>::max()),
        bestUpperobj(0),
        nzn(0),
        snzn(0),
        minLU(0),
        start_time(0),
        using_system_time(0),
        objSpeed(0) {
    lpsolvertime = shortestpathtime = 0;
    pivotType = NOCHANGE;
    linkId = -1;
  }
};

struct solverPara {
  LU_SOLVER solver;
  int maxIterationNum;
  int info;
  double objSpeedUpdateRat;
  int perIterationPrint;
  bool isSetpenaltyPrice;  // user set penalty price for fail demand
  double penaltyPriceForFailDemand;
  bool isSetDisturbed;
  int disturbedDigit;

  solverPara()
      : solver(LAPACK),
        maxIterationNum(1000000),
        info(0),
        objSpeedUpdateRat(0.4),
        perIterationPrint(100),
        isSetpenaltyPrice(false),
        isSetDisturbed(true),
        disturbedDigit(3) {}
};
}

#endif
