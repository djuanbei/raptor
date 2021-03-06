/**
 * @file   mcfgc.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Mon Apr  4 14:36:13 2016
 *
 * @brief  column generation method for multi-commodity flow problem
 *
 *
 */
#include <deque>
#include <vector>

#if (defined(__APPLE__) && defined(__MACH__))

#else
#include <omp.h>

#endif

#include <glpk.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include "System.h"
#include "config.h"
#include "graphalg.hpp"
#include "klu.h"
#include "sparse.h"
#include "util.h"

#ifdef CPLEX_SOLVER

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#define RC_EPS 1.0e-6

#endif

using namespace std;
using namespace raptor;
using namespace sparse;

using milliseconds = std::chrono::duration<double, std::milli>;

namespace raptor {

namespace mcmcf {

struct ENTER_VARIABLE {
  ENTER_BASE_TYPE type;
  int id;
  vector<int> path;
  ENTER_VARIABLE() : type(PATH_T), id(0) {}
};

struct Exit_base {
  EXIT_BASE_TYPE type;
  int id;
  Exit_base() : type(DEMAND_T), id(0) {}
};

struct KLUsolver {
  bool first;
  int cap;
  int *Ap;
  int *Ai;
  double *Ax;
  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common Common;
  int nonzeroNum;
  int dim;
  // vector<SparseMatrixElem> elements;
  KLUsolver() : first(true), Ap(NULL), Ai(NULL), Ax(NULL) {
#ifdef USING_KLU
    klu_defaults(&Common);
#endif
    nonzeroNum = 0;
    dim = 0;
    reScale(100);
  }
  ~KLUsolver();
  void setDim(int d) { dim = d; }
  void update(vector<SparseMatrixElem> &els, int n);
  bool solve(double *b);
  bool tsolve(double *b);
  void reScale(int s);
};

template <typename G, typename C, typename W = double>
class CG {
 public:
  struct Path {
    vector<int> path;
    int owner;
    int link;

    Path() : owner(-1), link(-1) {}
    Path(const vector<int> &p, const int o, const int l)
        : path(p), owner(o), link(l) {}
  };

 private:
  G graph;
  vector<int> linkMap;

  int origLink_num;

  W inf_weight;
  W inf;

  vector<Demand<C>> demands;

  vector<bool> orig_success;

  vector<W> orignal_weights;

  vector<C> orignal_caps;

  vector<W> update_weights;

  vector<C> update_caps;

  vector<C> edgeLeftBandwith;

  C TotalB;

  vector<Path> paths;  // all the paths save in this vector

  vector<Path> lastSP;

  vector<int> empty_paths;  // the location which  is delete path

  vector<double> min_commodity_cost;

  vector<int> saturate_links;  // contain the index of saturate links

  vector<int> saturate_link_ids;  // the id of  saturate link

  vector<int> un_saturate_links;

  vector<int> un_saturate_link_ids;  // the id of un saturate link

  vector<vector<int>>
      demand_secondary_path_locs;  // belong to fixed demands' second paths

  unordered_map<int, vector<int>>
      saturate_primary_path_locs;  // primary paths which
                                   // corross the saturate link

  vector<int> primary_path_loc;  // every demands have a primary path

  vector<int> saturate_link_path_loc;  //  the path corresponding saturate link

  Statistics_data sdata;

  vector<int> candidate_enters;

  KLUsolver klusolver;

  SparseSolver sparseSolver;


  double *rhs;
  double *X;
  double *Lambda;
  double *Mu;
  double * dual_solution; // it to value

  int *ipiv;
  double *b;

  double *SM;
  double *workS;
  int S_maxdim;

  double *x_K;
  double *x_S;
  double *x_N;

  double *lambda_K;
  double *lambda_S;
  double *lambda_N;

  C CZERO;
  /**
   * K: number of commodity
   * S: number of saturate link
   * N: number of unsaturate link
   */

  int K, S, N;
  double EPS;
  solverPara para;
  vector<pair<int, int>> saturateLinkAndMatrix;

  int thread_num;

  /**
   *
   *
   * @param S the size of saturate link
   */
  void allocateS(const int S) {
    if (S_maxdim >= S) {
      return;
    }
    if (NULL != SM) {
      delete[] SM;
      SM = NULL;
      delete[] workS;
    }
    S_maxdim = 1.3 * S + 1;
    SM = new double[S_maxdim * S_maxdim];
    workS = new double[S_maxdim * S_maxdim];
  }
  C success_obj() {
    C re = 0;
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      assert(pid < (int)paths.size() && !paths[pid].path.empty());
      if (paths[pid].path.back() < origLink_num) {
        re += x_K[i];
      }
    }

    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      int pid = saturate_link_path_loc[link];

      if (paths[pid].path.back() < origLink_num) {
        re += x_S[i];
      }
    }
    return re;
  }

  double computeOBJ() {
    double re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      double p_cost = path_cost(orignal_weights, paths[pid].path, 0.0);

      re += p_cost * x_K[i];
    }

    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      int pid = saturate_link_path_loc[link];

      double p_cost = path_cost(orignal_weights, paths[pid].path, 0.0);

      re += p_cost * x_S[i];
    }
    return re;
  }
  double computeObj() {
    double re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      if (paths[pid].path.back() < origLink_num) {
        double p_cost = path_cost(orignal_weights, paths[pid].path, 0.0);

        re += p_cost * x_K[i];
      }
    }

    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      int pid = saturate_link_path_loc[link];
      if (paths[pid].path.back() < origLink_num) {
        double p_cost = path_cost(orignal_weights, paths[pid].path, 0.0);

        re += p_cost * x_S[i];
      }
    }
    return re;
  }

  int nonzero(const int n, const double *data) {
    int re = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (fabs(data[i * n + j]) > EPS) {
          re++;
        }
      }
    }
    return re;
  }

  void deleteElem(vector<int> &data, int key) {
    vector<int>::iterator it = find(data.begin(), data.end(), key);
    data.erase(it);
  }
  /**
   * row primary
   *
   */

  void getSparseS(vector<SparseMatrixElem> &elements) {
    /**
     * -D
     */
    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      int pid = saturate_link_path_loc[link];
      const vector<int> &path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        int j = saturate_link_ids[*lit];
        if (j > -1) {
          SparseMatrixElem elem;
          elem.column = i;
          elem.row = j;
          elem.value = -1.0;
          elements.push_back(elem);
        }
      }
    }

    /**
     *  CB
     *
     */
    for (unordered_map<int, vector<int>>::const_iterator it =
             saturate_primary_path_locs.begin();
         it != saturate_primary_path_locs.end(); it++) {
      int link = it->first;
      int i = saturate_link_ids[link];  
      assert(i > -1);
      for (vector<int>::const_iterator pit = it->second.begin();
           pit != it->second.end(); pit++) {
        int oindex = paths[*pit].owner;
        for (vector<int>::const_iterator cit =
                 demand_secondary_path_locs[oindex].begin();
             cit != demand_secondary_path_locs[oindex].end(); cit++) {
          int slink = paths[*cit].link;
          int j = saturate_link_ids[slink]; 
          assert(j > -1);

          SparseMatrixElem elem;
          elem.column = j;
          elem.row = i;
          elem.value = 1.0;
          elements.push_back(elem);
        }
      }
    }
  }

  void computeS() {
    if (0 == S) {
      return;
    }
    double t = 0;
    if (KLU == para.solver) {
      vector<SparseMatrixElem> elements;
      getSparseS(elements);

      callTime(t, klusolver.update(elements, S));

      sdata.nzn = klusolver.nonzeroNum;

    } else if (LAPACK == para.solver) {
      /**
       * column primary
       *
       */
      /**
       * SM= CB - D
       *
       */
      sdata.nzn = 0;
      if (S > S_maxdim) {
        allocateS(S);
      }
      fill(SM, SM + S * S, 0.0);
      /**
       * -D
       */
      for (int i = 0; i < S; i++) {
        int link = saturate_links[i];
        int pid = saturate_link_path_loc[link];
        const vector<int> &path = paths[pid].path;
        for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
             lit++) {
          int j = saturate_link_ids[*lit];
          if (j > -1) {
            SM[j * S + i] = -1.0;
            sdata.nzn++;
          }
        }
      }

      /**
       *  CB
       *
       */
      for (unordered_map<int, vector<int>>::const_iterator it =
               saturate_primary_path_locs.begin();
           it != saturate_primary_path_locs.end(); it++) {
        int link = it->first;
        int i = saturate_link_ids[link]; 
        assert(i > -1);
        for (vector<int>::const_iterator pit = it->second.begin();
             pit != it->second.end(); pit++) {
          int oindex = paths[*pit].owner;
          for (vector<int>::const_iterator cit =
                   demand_secondary_path_locs[oindex].begin();
               cit != demand_secondary_path_locs[oindex].end(); cit++) {
            int slink = paths[*cit].link;
            int j = saturate_link_ids[slink];
            SM[i * S + j] += 1.0;
            sdata.nzn++;
          }
        }
      }
    } else if (SPARSE == para.solver) {
      vector<SparseMatrixElem> elements;
      getSparseS(elements);
      callTime(t, sparseSolver.update(elements));
      sdata.nzn = sparseSolver.getNonzero();
    }
    sdata.lpsolvertime += t;
  }

  void transposeS() {
    double temp = 0.0;
    for (int i = 0; i < S - 1; i++) {
      for (int j = i + 1; j < S; j++) {
        temp = SM[i * S + j];
        SM[i * S + j] = SM[j * S + i];
        SM[j * S + i] = temp;
      }
    }
  }

  /**
   * Ax=rhs
   *
   * [  I_{K*K}   B            0       ]  [ x_K ] = [ rhs_K ]
   * [  C         D            0       ]  [ x_S ] = [ rhs_S ]
   * [  H         F            I_{N*N} ]  [ x_N ] = [ rhs_N ]
   *
   *
   * x_S=(B rhs_K- rhs_S) /( CB-D )
   *
   * [  I_{K*K}   B           0       ]  [ x_K ] = [ rhs_K ]  bandwidth
   *envery
   *demand
   * [  C         D           0       ]  [ x_S ] = [ rhs_S ]  capacity of
   *saturate
   *links
   * [  H         F           I_{N*N} ]  [ x_N ] = [ rhs_N ]  capacity of
   *other
   *links
   *rhs  =  [rhs_K  rhs_S   rhsa_N]
   *x_N  =  ( C rhs_K -rhs_S )/(CB-D)=( C rhs_K -rhs_S )/(SM)
   **/
  void computeRHS() {
    x_K = X;
    x_S = X + K;
    x_N = X + K + S;

    int nrhs = 1;
    int lda = S;

    int ldb = S;
    int info = 0;

    if (S > 0) {
      fill(b, b + S, 0.0);
      /**
       * -rhs_S
       *
       */
      int i = 0;
      for (i = 0; i + 3 < S; i += 4) {
        b[i] = -rhs[K + saturate_links[i]];
        b[i + 1] = -rhs[K + saturate_links[i + 1]];
        b[i + 2] = -rhs[K + saturate_links[i + 2]];
        b[i + 3] = -rhs[K + saturate_links[i + 3]];
      }
      for (; i < S; i++) {
        b[i] = -rhs[K + saturate_links[i]];
      }
      /**
       * C rhs_K
       *
       */

      for (int i = 0; i < S; i++) {
        int link = saturate_links[i];

        const vector<int> &pps = saturate_primary_path_locs[link];
        for (vector<int>::const_iterator it = pps.begin(); it != pps.end();
             it++) {
          b[i] += rhs[paths[*it].owner];
        }
      }

      double t = 0;

      if (KLU == para.solver) {
        callTime(t, klusolver.solve(b));

      } else if (LAPACK == para.solver) {
        copy(SM, SM + S * S, workS);
        callTime(t, dgesv_(&S, &nrhs, (double *)workS, &lda, ipiv, (double *)b,
                           &ldb, &info));
        sdata.snzn = nonzero(S, workS);

        if (info > 0) {
          assert(false);
          printf(
              "The diagonal element of the "
              "triangular factor of "
              "Lambda,\n");
          printf(
              "U(%i,%i) is zero, so that Lambda is "
              "singular;\n",
              info, info);
          printf(
              "the solution could not be "
              "computed.\n");
          exit(1);
        }
      } else if (SPARSE == para.solver) {
        callTime(t, sparseSolver.incSolver(x_S, b));
      } else {
        assert(false);
      }
      sdata.lpsolvertime += t;
      memcpy(x_S, b, S * sizeof(double));
#ifdef DEBUG
      for (int i = 0; i < S; i++) {
        int link = saturate_links[i];
        int spid = saturate_link_path_loc[link];
        assert(b[i] <= rhs[paths[spid].owner]+EPS);
      }
#endif
    }

    /**
     * x_K=rhs_K-B x_S
     *
     */
    memcpy(x_K, rhs, K * sizeof(double));

    for (int i = 0; i < K; i++) {
      const vector<int> &pathindices = demand_secondary_path_locs[i];
      for (vector<int>::const_iterator it = pathindices.begin();
           it != pathindices.end(); it++) {
        x_K[i] -= x_S[saturate_link_ids[paths[*it].link]];
      }
    }
    /**
     * x_N=rhs_N- H x_K-F x_S
     *
     */
    for (int i = 0; i < N; i++) {
      int link = un_saturate_links[i];
      x_N[i] = rhs[K + link];
    }
    computeLeftN(X);

#ifdef DEBUG
    for (int i = 0; i < K + N + S; i++) {
      assert(X[i] >= -EPS);
    }
#endif
  }

  /**
   *
   * @brief choose candiate exit base
   *
   * @return
   */
  Exit_base computeExitVariable() {
    /**
     *  choose entering base as i=argmin_{ lambda[ i]>0} x[ i ]/lambda[ i ]
     *
     */

    int reK = -1;
    double minK = numeric_limits<double>::max();
    double temp;
    int i = 0;
    
    int num=K+S+N;

    while (i < num) {
      if (Lambda[i] > EPS) {
        temp = X[i] / Lambda[i];
        if (temp < minK) {
          reK = i;
          minK = temp;
        }
      }
      i++;
    }
    assert(reK>-1);
    if (minK <= EPS / 100.0) {
        sdata.empty_iterator_num++;
    }
    Exit_base re;
    if(reK<K){
      re.id = reK;
      re.type = DEMAND_T;
    }else if(reK<K+S){
      re.id = saturate_links[reK-K];
      re.type = STATUS_LINK;
    }else{
      re.id = un_saturate_links[reK-K-S];
      re.type = OTHER_LINK;
    }
    return re;

  }

 public:
  CG(const G &g, const vector<W> &ws, const vector<C> &caps,
     const vector<Demand<C>> &ds)
      : demands(ds), thread_num(1) {
    orignal_weights.resize(ws.size());
    orignal_caps.resize(ws.size());
    vector<int> srcs, snks;
    linkMap.resize(ws.size());
    rhs = NULL;
    /**
     *  Reorder links' id to improve
     * graph get adjacent nodes speed
     *
     */

    for (size_t v = 0; v < g.getVertex_num(); v++) {
      int degree = g.getOutDegree(v);
      for (int i = 0; i < degree; i++) {
        int link = g.getAdj(v, i);
        int id = srcs.size();
        linkMap[id] = link;
        orignal_caps[id] = caps[link];
        orignal_weights[id] = ws[link];
        int snk;
        g.findRhs(link, v, snk);
        srcs.push_back(v);
        snks.push_back(snk);
      }
    }
    graph.initial(srcs, snks);

    double max_w = 0;
    origLink_num = graph.getLink_num();
    for (int i = 0; i < origLink_num; i++) {
      if (max_w < orignal_weights[i]) {
        max_w = orignal_weights[i];
      }
    }
    if (para.isSetpenaltyPrice) {
      inf_weight = para.penaltyPriceForFailDemand;

    } else {
      inf_weight = max_w * graph.getVertex_num() + 1;
    }

    inf = inf_weight;

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < 1000; i++) {
      thread_num = omp_get_num_threads();
    }

#endif  // _OPENMP

    saturate_link_path_loc.resize(origLink_num + 1);

    CZERO = ((C)1e-6);
    EPS = ((double)1e-4);

    K = demands.size();
    demand_secondary_path_locs.resize(K);

    S = 0;
    N = origLink_num + 1;

    b = new double[N];
    X = new double[K + N];
    Lambda = new double[K + N];
    Mu = new double[N];
    dual_solution =new double[N];
    fill(Mu, Mu + N, 0.0);

    fill(X, X + N, 0.0);
    ipiv = new int[N];

    for (int i = 0; i < K; i++) {
      X[i] = demands[i].bandwidth;
    }

    x_K = X;
    lambda_K = Lambda;

    SM = NULL;
    workS = NULL;
    S_maxdim = 0;

    update_caps = orignal_caps;

    C minNonzero = numeric_limits<C>::max();
    TotalB = 0;
    for (typename vector<Demand<C>>::const_iterator it = demands.begin();
         it != demands.end(); it++) {
      if (it->bandwidth > CZERO && it->bandwidth < minNonzero) {
        minNonzero = it->bandwidth;
      }
      TotalB += it->bandwidth;
    }

    // if (para.isSetDisturbed) {
    //   double domain = pow(10, para.disturbedDigit) + 1;
    //   double multiRat = 1.0 / domain;

    //   if (minNonzero > 1) {
    //     for (size_t i = 0; i < update_caps.size(); i++) {
    //       update_caps[i] += multiRat * (rand() / (RAND_MAX + 0.0));
    //     }

    //   } else {
    //     for (size_t i = 0; i < update_caps.size(); i++) {
    //       update_caps[i] += multiRat * minNonzero * (rand() / (RAND_MAX + 0.0));
    //     }
    //   }
    // }
    /**
     *add a  dummy link which can setup all demands
     *
     */

    orignal_weights.push_back(inf_weight);
    update_caps.push_back(TotalB + 1);
    inf_weight = getInf(0.0);
#ifdef STATIC_TABLE
    int did;
    getData(graph.getVertex_num(), inf_weight, did);
    getData(graph.getVertex_num(), inf_weight, did, false);
#endif
  }
  ~CG() {
    if (NULL != Lambda) {
      delete[] Lambda;
      Lambda = NULL;
    }
    if (NULL != rhs) {
      delete[] rhs;
      rhs = NULL;
    }
    if (NULL != X) {
      delete[] X;
      X = NULL;
    }
    if (NULL != Mu) {
      delete[] Mu;
      Mu = NULL;
    }
    if(NULL!=dual_solution){
      delete[] dual_solution;
      dual_solution=NULL;
      
    }
    if (NULL != ipiv) {
      delete[] ipiv;
      ipiv = NULL;
    }

    if (NULL != b) {
      delete[] b;
      b = NULL;
    }

    if (NULL != SM) {
      delete[] SM;
      SM = NULL;
      delete[] workS;
      workS = NULL;
    }
  }
  void setPara(solverPara &p) { para = p; }

  void setInfo(const int level) { para.info = level; }

  void setLUSOLVER(LU_SOLVER s) { para.solver = s; }

  void setStatusLink(const int link, const int pid) {
    paths[pid].link = link;
    saturate_link_path_loc[link] = pid;
  }

  void addPrimarySaturateLink(const int commodityID) {
    int pid = primary_path_loc[commodityID];
    const vector<int> &new_path = paths[pid].path;
    for (vector<int>::const_iterator it = new_path.begin();
         it != new_path.end(); it++) {
      if (saturate_link_ids[*it] > -1) {
        saturate_primary_path_locs[*it].push_back(pid);
      }
    }
  }

  void deletePrimarySatuateLinks(const int commodityID) {
    int primary_pid = primary_path_loc[commodityID];
    const vector<int> &orig_path = paths[primary_pid].path;
    for (vector<int>::const_iterator it = orig_path.begin();
         it != orig_path.end(); it++) {
      if (saturate_link_ids[*it] > -1) {
        deleteElem(saturate_primary_path_locs[*it], primary_pid);

      }
    }
  }

  void addStatusLink(const int link) {
    saturate_primary_path_locs[link].clear();
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      if (find(paths[pid].path.begin(), paths[pid].path.end(), link) !=
          paths[pid].path.end()) {
        saturate_primary_path_locs[link].push_back(pid);
      }
    }
  }
  void deleteSaturateLink(const int link) {
    vector<int>::iterator it =
        find(saturate_links.begin(), saturate_links.end(), link);
#if DEBUG
    assert(*it == link && it != saturate_links.end());
#endif

    saturate_links.erase(it);
    if (SPARSE == para.solver) {    
      for (int i = saturate_link_ids[link]; i + 1 < S; i++) {
        x_S[i] = x_S[i + 1];
        Mu[i] = Mu[i + 1];
      }
    }

    saturate_link_ids[link] = -1;

    int spid = saturate_link_path_loc[link];
    paths[spid].link = -1;

    saturate_primary_path_locs.erase(link);
  }

  void changeSaturateLink(const int link1, const int link2) {
    vector<int>::iterator it =
        find(saturate_links.begin(), saturate_links.end(), link1);
    *it = link2;

    saturate_link_ids[link2] = saturate_link_ids[link1];

    saturate_link_ids[link1] = -1;

    int spid = saturate_link_path_loc[link1];
    paths[spid].link = link2;

    saturate_primary_path_locs.erase(link1);

    saturate_link_path_loc[link2] = spid;

    addStatusLink(link2);
  }

  double getOrigCost(const vector<int> &path) const {
    return path_cost(orignal_weights, path, 0.0);
  }

  void computIndexofLinks() {

    for (size_t i = 0; i < saturate_links.size(); i++) {
      saturate_link_ids[saturate_links[i]] = i;
    }

    int i = 0;

    for (int j = 0; j < S + N; j++) {
      if (saturate_link_ids[j] < 0) {
        un_saturate_link_ids[j] = i;
        un_saturate_links[i++] = j;
      } else {
        un_saturate_link_ids[j] = -1;
      }
    }
  }

#ifdef CPLEX_SOLVER

  void report3(IloCplex &masterSolver, int start, double penalty,
               IloNumVarArray usedVariable) {
    double OBJ = masterSolver.getObjValue();
    for (int i = 0; i < start; i++) {
      OBJ -= masterSolver.getValue(usedVariable[i]) * penalty;
    }
    cout << "Minimum cost: " << OBJ << endl;

    double succ = 0;

    for (IloInt j = start; j < usedVariable.getSize(); j++) {
      succ += masterSolver.getValue(usedVariable[j]);
    }
    cout << "Maximum success bw: " << succ << endl;
  }

  bool CPLEX_solve() {
    IloEnv env;
    try {
      IloInt i, j;
      IloNum penalty = inf;

      IloNumArray bws(env);

      IloNumArray dbws(env);

      for (size_t i = 0; i < orignal_caps.size(); i++) {
        bws.add(orignal_caps[i]);
      }

      for (size_t k = 0; k < demands.size(); k++) {
        dbws.add(demands[k].bandwidth);
      }
      int commodityNum = demands.size();
      bws.add(TotalB + 1);

      /// MASTER-OPTIMIZATION PROBLEM ///

      IloModel mcfModel(env);
      IloObjective totalPrice = IloAdd(mcfModel, IloMinimize(env));

      IloRangeArray demandCons =
          IloAdd(mcfModel, IloRangeArray(env, dbws, IloInfinity));
      IloRangeArray linkCons = IloAdd(mcfModel, IloRangeArray(env, 0, bws));

      IloNumVarArray usedVariable(env);

      IloInt linkNum = bws.getSize() - 1;

      vector<vector<int>> lastSP(commodityNum);

      for (j = 0; j < commodityNum; j++) {
        usedVariable.add(IloNumVar(totalPrice(penalty) + demandCons[j](1) +
                                   linkCons[linkNum](1)));
      }

      IloCplex masterSolver(mcfModel);

      /// COLUMN-GENERATION PROCEDURE ///
      vector<double> price(linkNum);

      IloNumArray newPathVarible(env, linkNum);

      /// COLUMN-GENERATION PROCEDURE ///
      bool state = true;
      for (; state;) {
        /// OPTIMIZE OVER CURRENT PATTERNS ///
        state = false;
        masterSolver.solve();

        for (i = 0; i < linkNum; i++) {
          price[i] = orignal_weights[i] - masterSolver.getDual(linkCons[i]);
        }
        vector<int> path;
        for (int i = 0; i < commodityNum; i++) {
          /// FIND AND ADD A NEW PATTERN ///
          if (bidijkstra_shortest_path(graph, price, demands[i].src,
                                       demands[i].snk, path, inf_weight)) {
            if (lastSP[i].empty() || (path_cost(price, lastSP[i], 0.0) >
                                      path_cost(price, path, 0.0))) {
              lastSP[i] = path;
              state = true;
              for (int i = 0; i < linkNum; i++) {
                newPathVarible[i] = 0;
              }
              for (vector<int>::iterator it = path.begin(); it != path.end();
                   it++) {
                newPathVarible[*it] = 1;
              }
              double unitPrice = path_cost(orignal_weights, path, 0.0);
              usedVariable.add(IloNumVar(totalPrice(unitPrice) +
                                         demandCons[i](1) +
                                         linkCons(newPathVarible)));
            }
          }
        }
      }
      cout << "Solution status: " << masterSolver.getStatus() << endl;
      double pp = penalty;
      report3(masterSolver, commodityNum, pp, usedVariable);
    } catch (IloException &ex) {
      cerr << "Error: " << ex << endl;
    } catch (...) {
      cerr << "Error" << endl;
    }
    env.end();

    return true;
  }

#endif

  bool GlkpSolve() {
    glp_term_out(GLP_OFF);
    // Creates the GLPK problem instance.
    glp_prob *lp = glp_create_prob();
    glp_set_prob_name(lp, "Multi-commodity flow problem");
    // The objective function will be minimized.
    glp_set_obj_dir(lp, GLP_MIN);
    // Adds the constraints (rows)
    int commodityNum = demands.size();
    int linkNum = orignal_caps.size();

    glp_add_rows(lp, commodityNum + linkNum + 1);
    for (int i = 1; i <= commodityNum; i++) {
      glp_set_row_bnds(lp, i, GLP_FX, demands[i - 1].bandwidth, 0);
    }
    for (int i = 1; i <= linkNum; i++) {
      glp_set_row_bnds(lp, i + commodityNum, GLP_UP, 0, update_caps[i - 1]);
    }
    glp_set_row_bnds(lp, commodityNum + linkNum + 1, GLP_DB, 0, TotalB + 1);

    // Adds the variables (columns) and also set the coefficients of the
    // objective function.
    glp_add_cols(lp, commodityNum);

    for (int i = 1; i <= commodityNum; i++) {
      glp_set_col_bnds(lp, i, GLP_LO, 0, 0);  // x[i] >= 0
      glp_set_obj_coef(lp, i, inf);           // c[i] = inf
    }
    int m = commodityNum + linkNum + 1;
    const int size = 2 * commodityNum + 2;

    int A_row[size], A_column[size];
    double A_value[size];
    int index = 1;
    for (int j = 1; j <= commodityNum; j++) {
      // Adds this intial column to A
      A_column[index] = j;   // column j
      A_row[index] = j;      // line j
      A_value[index] = 1.0;  // a[j,j] = 1
      index++;

      A_column[index] = j;  // column j
      A_row[index] =
          commodityNum + linkNum + 1;  // line commodityNum + linkNum + 1
      A_value[index] = 1.0;            // a[commodityNum + linkNum + 1,j] = 1
      index++;
    }

    // Loads the A matrix.
    glp_load_matrix(lp, index - 1, A_row, A_column, A_value);

    int n = commodityNum;
    vector<vector<int>> lastSP(commodityNum);
    vector<double> price(linkNum);
    //
    // PROBLEM SOLVING ITERATIONS
    //
    int columnIndex[3 + linkNum];
    double columnValue[3 + linkNum];
    bool stop = false;


    do {
      stop = true;
      // Solve the current iteration's linear program using SIMPLEX.
      int statusCode = glp_simplex(lp, NULL);
      // The SIMPLEX return code should be 0; means no errors.
      assert(statusCode == 0);

      // Gets the dual variables

      for (int i = 0; i < linkNum; i++) {
        price[i] =
            orignal_weights[i] - glp_get_row_dual(lp, i + commodityNum + 1);
      }

      vector<int> path;
      for (int i = 0; i < commodityNum; i++) {
        /// FIND AND ADD A NEW PATTERN ///
        if (bidijkstra_shortest_path(graph, price, demands[i].src,
                                     demands[i].snk, path, inf_weight)) {
          if (lastSP[i].empty() || (path_cost(price, lastSP[i], 0.0) >
                                    path_cost(price, path, 0.0))) {

            lastSP[i] = path;
            stop = false;
            columnIndex[1] = i + 1;
            columnValue[1] = 1.0;

            for (size_t j = 0; j < path.size(); j++) {
              columnIndex[j + 2] = commodityNum + 1 + path[j];
              columnValue[j + 2] = 1.0;
            }
            glp_add_cols(lp, 1);
            int j = glp_get_num_cols(lp);
            double unitPrice = path_cost(orignal_weights, path, 0.0);
            glp_set_col_bnds(lp, j, GLP_LO, 0, 0);
            glp_set_obj_coef(lp, j, unitPrice);

            glp_set_mat_col(lp, j, 1 + path.size(), columnIndex, columnValue);

            n = n + 1;
          }
        }
      }

      assert(n == glp_get_num_cols(lp));
    } while (not stop);

    double OBJ = glp_get_obj_val(lp);

    for (int i = 1; i <= commodityNum; i++) {
      OBJ -= glp_get_col_prim(lp, i) * inf;
    }
    cout << fixed;
    cout << "Minimum cost: " << OBJ << endl;

    double succ = 0;

    for (int j = commodityNum + 1; j < n + 1; j++) {
      succ += glp_get_col_prim(lp, j);
    }
    cout << "Maximum success bw: " << succ << endl;
    cout<<"variable num: "<<n-commodityNum<<endl;
    glp_delete_prob(lp);
    return true;
  }

  bool solve() {
    initial_solution();
    sdata.start_time = systemTime();
    GlkpSolve();
    cout << "total take " << systemTime() - sdata.start_time << endl;

#ifdef CPLEX_SOLVER
    CPLEX_solve();
    cout << "total take " << systemTime() - sdata.start_time << endl;

#endif
    sdata.start_time = systemTime();
    bool re = iteration();

    if (para.info > 0) {
      printResult();
    }
    return re;
  }
  void writeKsptoCNF(const int k, const char *outfile) {
    initial_solution();
    vector<C> temp_cap(orignal_caps);
    vector<vector<vector<int>>> paths(K);
    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      vector<W> ws = orignal_weights;
      for (int i = 0; i < origLink_num; i++) {
        if (orignal_caps[i] < bw) {
          ws[i] = 2 * inf_weight;
        }
      }

      inc_ksp::yen_ksp<G, vector<W>, W> ksp_g(graph, orignal_weights,
                                              inf_weight);
      inc_ksp::yen_next_path<G, vector<W>, W> ksp = ksp_g.next_path(src, snk);
      vector<int> path;
      for (int j = 0; j < k && ksp.next_path(path); j++) {
        paths[i].push_back(path);
      }
    }
  }

  /**
   * Use some heustic  to choose good initial flow for every commodity.
   *Set a initial routing path for every demand, the only constraint the
   * initial solution must satisfies is that link capacity.
   *
   */

  void initial_solution() {
    paths.resize(K);
    primary_path_loc.resize(K, 0);
    orig_success.resize(K, true);

/**
 *  check whether there is path from  demand source to demand target
 *
 */

#pragma omp parallel for

    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      vector<int> path;
      if (!bidijkstra_shortest_path(graph, orignal_weights, src, snk, path,
                                    inf)) {
        orig_success[i] = false;
      }
    }

    vector<C> temp_cap(orignal_caps);
    vector<bool> success(K, false);

    vector<pair<C, int>> sortDemands;
    for (int i = 0; i < K; i++) {
      pair<C, int> temp = make_pair(demands[i].bandwidth, i);
      sortDemands.push_back(temp);
    }

    sort(sortDemands.rbegin(), sortDemands.rend());

    /**
     * Direcltly setup demand one by one
     *
     */

    for (int k = 0; k < K; k++) {
      int i = sortDemands[k].second;

      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      vector<W> ws(orignal_weights.size(), 1);

      for (int j = 0; j < origLink_num; j++) {
        if (temp_cap[j] < bw) {
          ws[j] = graph.getVertex_num();
        } else {
          ws[j] = bw / temp_cap[j];
        }
      }

      vector<int> path;
      if (bidijkstra_shortest_path(graph, ws, src, snk, path,
                                   graph.getVertex_num())) {
        success[i] = true;

        for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
          temp_cap[*it] -= bw;
        }

        paths[i].path.swap(path);
        paths[i].owner = i;
        primary_path_loc[i] = i;
      }
    }

    for (int i = 0; i < K; i++) {
      if (!success[i]) {
        vector<int> path;
        path.push_back(origLink_num);
        paths[i].path.swap(path);
        paths[i].owner = i;
        primary_path_loc[i] = i;
        success[i] = true;
      }
    }

    update_weights = orignal_weights;

    min_commodity_cost.resize(K, 2 * inf_weight);
#pragma omp parallel for
    for (int i = 0; i < K; i++) {
      vector<int> path;
      int src = demands[i].src;
      int snk = demands[i].snk;
      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                                   inf_weight)) {
        min_commodity_cost[i] = path_cost(update_weights, path, (W)0.0);
      }
    }
    saturate_link_ids.resize(N, -1);

    for (int i = 0; i < N; i++) {
      un_saturate_links.push_back(i);
      un_saturate_link_ids.push_back(i);
    }
    rhs = new double[K + origLink_num + 1];
    // rhs.resize(K + origLink_num + 1, (C)0.0);
    for (int i = 0; i < K; i++) {
      rhs[i] = demands[i].bandwidth;
    }

    for (size_t i = 0; i < update_caps.size(); i++) {
      rhs[i + K] = update_caps[i];
    }

    // dual_solution.resize(N + 1, 0);
  }

  bool iteration() {
    vector<ENTER_VARIABLE>  enterVariables;
    
    computeRHS();
    while (true) {
      sdata.iterator_num++;


      if (para.info > 2) {
        double OBJ = computeOBJ();
        sdata.totalStaturateLink += S;
        sdata.totalNonzero += sdata.nzn;
        if (para.info > 1) {
          saturateLinkAndMatrix.push_back(make_pair(S, sdata.nzn));
        }

        if (sdata.bestUpperobj < OBJ - sdata.estimee_opt_diff) {
          sdata.bestUpperobj = OBJ - sdata.estimee_opt_diff;
        }
        if (sdata.iterator_num % para.perIterationPrint == 0) {
          sdata.using_system_time = systemTime() - sdata.start_time;
          C sobj = success_obj();

          std::cout << fixed;
          std::cout << "============================================ " << endl;
          std::cout << "iteration: " << sdata.iterator_num << endl;

          std::cout << "using time(s) :" << sdata.using_system_time
                    << std::endl;
          std::cout << "empty iteration nonzero_bata: "
                    << sdata.empty_iterator_num << endl;
          std::cout << "saturate link nonzero_bata: " << S << endl;
          std::cout << "objvalue: " << OBJ << endl;
          std::cout << "success fractional bw: " << sobj
                    << ", success rat: " << sobj / (TotalB + 0.01) << std::endl;
          std::cout << "The obj gap from opt less or equal than: "
                    << sdata.estimee_opt_diff << "("
                    << 100 * sdata.estimee_opt_diff / sdata.bestUpperobj << "%)"
                    << std::endl;

          std::cout << "The number of nonzeon element in matrix (CB-D) : "
                    << sdata.nzn << std::endl;

          std::cout << "Last entering: " << sdata.enter
                    << ", last exit: " << sdata.exit << std::endl;
        }
      }
      if (sdata.iterator_num > para.maxIterationNum) {
        return false;
      }

      // update_edge_left_bandwith();
      /**
       *  entering variable choose
       *
       */
      double t = 0;
      callTime(t,  chooseEnteringVariable(enterVariables));
      sdata.shortestpathtime += t;

      if (enterVariables.empty()) {
        return true;
      }

      Exit_base exit_base;
      /**
       *  exit base  choose
       *
       */
      if (PATH_T == enterVariables[0].type) {
        exit_base = getExitBasebyPath(enterVariables[0]);
      } else {
        exit_base = getExitBasebyStatusLink(enterVariables[0]);
      }

      pivot(enterVariables[0], exit_base);

      S = saturate_links.size();
      assert(empty_paths.size() + K + S == paths.size());
      N = origLink_num + 1 - S;
      computIndexofLinks();

      computeS();

#pragma omp sections 
      {
#pragma omp section
        {

          update_edge_cost();
        }
        
#pragma omp section
        {
          /**
           * column primary
           *
           */
          if (LAPACK == para.solver) {
            transposeS();
          }
          computeRHS();
        }
      }
    }
  }
  /**
   *
   * choose a comodity which the current best path has biggest diff from
   *old
   *solution
   *
   * @return
   */
  void chooseEnteringVariable(vector<ENTER_VARIABLE> & enterVariables) {
    enterVariables.clear();
    ENTER_VARIABLE enter_variable;
    enter_variable.type = LINK_T;
    enter_variable.id = -1;
    /**
     *  check whether saturate link dual value is negative
     *
     */
    double max_diff = -EPS;

    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      if (dual_solution[i] < max_diff) {
        max_diff = dual_solution[i];
        enter_variable.id = link;
      }
    }

    /**
     * If there is a dual value of  saturate link is negative then this link
     * is
     * a entering variable
     */

    if (enter_variable.id >= 0) {
      enterVariables.push_back(enter_variable);
      return;
    }

    /**
     * fast way to choose entering variable as a path
     */
    for (size_t i = 0; i < candidate_enters.size(); i++) {
      int id = candidate_enters[i];
      int src = demands[id].src;
      int snk = demands[id].snk;
      vector<int> path;

      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[id]].path, (W)0.0);

      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                                   inf_weight)) {
        W new_cost = path_cost(update_weights, path, (W)0.0);
        W diff = old_cost - new_cost;

        if (diff > 1000 * EPS && diff > sdata.objSpeed){

          candidate_enters.erase(candidate_enters.begin(),
                                candidate_enters.begin() + i + 1);
          enter_variable.id = id;
          enter_variable.type = PATH_T;
          enter_variable.path.swap(path);

          sdata.objSpeed = (para.objSpeedUpdateRat) * sdata.objSpeed +
                           (1 - (para.objSpeedUpdateRat)) * diff;
          enterVariables.push_back(enter_variable);
          return;
        }
      }
    }

    candidate_enters.clear();

    sdata.estimee_opt_diff = 0;
    vector<double> opt_gap(thread_num, 0);
    vector<double> max_diffs(thread_num, EPS);

    vector<ENTER_VARIABLE> enter_variables(thread_num);

    int chid = -1;
    max_diff = EPS;
    for (int i = 0; i < thread_num; i++) {
      enter_variables[i].type = PATH_T;
      enter_variables[i].id = -1;
    }
    // cout<<"lastSP size: "<<lastSP.size()<<endl;
#pragma omp parallel for
    for(size_t k=0; k< lastSP.size(); k++){

#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;

#endif  // (_OPENMP)
      int i=lastSP[k].owner;
      if (!orig_success[i]) {
        continue;
      }
      int src = demands[i].src;
      int snk = demands[i].snk;

      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]].path, (W)0.0);
      W new_cost = path_cost(update_weights, lastSP[k].path, (W)0.0);
      W temp_diff = (old_cost - new_cost);
      if(temp_diff>EPS){
        if(temp_diff>max_diffs[tid]){
          max_diffs[tid]=temp_diff;
          enter_variables[tid].id = i;
          enter_variables[tid].path=lastSP[k].path;
        }
      }
    }
    
    chid = 0;
    max_diff = max_diffs[0];

    for (int i = 1; i < thread_num; i++) {
      if (max_diffs[i] > max_diff) {
        max_diff = max_diffs[i];
        chid = i;
      }
    }
    
    if(enter_variables[chid].id>-1){
      enterVariables.push_back(enter_variables[chid]);
      return ;
    }



    vector<vector<int>> candidate(thread_num);
    vector<vector<int>> path(thread_num);
#pragma omp parallel for
    for (int i = 0; i < K; i++) {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;

#endif  // (_OPENMP)
      if (!orig_success[i]) {
        continue;
      }
      int src = demands[i].src;
      int snk = demands[i].snk;

      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]].path, (W)0.0);
      /**
       * If the possible biggest reduce objective value is less than exist
       * found
       * one then continue
       *
       */

      if ((old_cost - min_commodity_cost[i]) * demands[i].bandwidth <
          max_diffs[tid]) {
        continue;
      }

      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path[tid],
                                   inf_weight)) {
        W new_cost = path_cost(update_weights, path[tid], (W)0.0);

        if (new_cost > inf) {
          path[tid].clear();
          path[tid].push_back(origLink_num);
          new_cost = inf;
        }

        W temp_diff = (old_cost - new_cost);
        // if (lastSP[i] != path[tid]) {
        //   lastSP[i] = path[tid];
        // }

        if (temp_diff > EPS) {
          Path temp;
          temp.owner=i;
          temp.path=path[tid];
          lastSP.push_back(temp);
          if (temp_diff > max_diffs[tid]) {
            max_diffs[tid] = temp_diff;
            opt_gap[tid] += temp_diff * demands[i].bandwidth;

            if (temp_diff > 10000 * EPS) {
              candidate[tid].push_back(i);
            }

            enter_variables[tid].id = i;
            enter_variables[tid].path.swap(path[tid]);
          }
        }
      }
    }

    chid = 0;
    max_diff = max_diffs[0];
    candidate_enters = candidate[0];
    for (int i = 1; i < thread_num; i++) {
      candidate_enters.insert(candidate_enters.end(), candidate[i].begin(),
                             candidate[i].end());
      if (max_diffs[i] > max_diff) {
        max_diff = max_diffs[i];
        chid = i;
      }
    }

    for (int i = 0; i < thread_num; i++) {
      sdata.estimee_opt_diff += opt_gap[i];
    }

    sdata.objSpeed = (para.objSpeedUpdateRat) * sdata.objSpeed +
                     (1 - (para.objSpeedUpdateRat)) * max_diffs[chid];
    if(enter_variables[chid].id>-1){
      enterVariables.push_back(enter_variables[chid]);
    }
  }

  /**
   * data_N -=  -H data_K-F data_S
   *
   */
  void computeLeftN(double *data) {
    double *data_K = data;
    double *data_S = data + K;
    double *data_N = data + K + S;

    for (int pid = 0; pid < K; pid++) {
      if (data_K[pid] != 0) {
        int ppid = primary_path_loc[pid];
        for (vector<int>::const_iterator lit = paths[ppid].path.begin();
             lit != paths[ppid].path.end(); lit++) {
          int i = un_saturate_link_ids[*lit];

          if (i > -1) {
            data_N[i] -= data_K[pid];
          }
        }
      }
    }

    for (int lid = 0; lid < S; lid++) {
      if (data_S[lid] != 0) {
        int ppid = saturate_link_path_loc[saturate_links[lid]];
        for (vector<int>::const_iterator lit = paths[ppid].path.begin();
             lit != paths[ppid].path.end(); lit++) {
          int i = un_saturate_link_ids[*lit];

          if (i > -1) {
            data_N[i] -= data_S[lid];
          }
        }
      }
    }
  }

  /**
   * A lambda=beta
   *
   * [  I_{K*K}   B            0       ]  [ lambda_K ] = [ beta_K ]
   * [  C         D            0       ]  [ lambda_S ] = [ beta_S ]
   * [  H         F            I_{N*N} ]  [ lambda_N ] = [ beta_N ]

   * lambda_S=(B beta_K- beta_S) /( CB-D )
   *beta  =  [beta_K  beta_S   beta_N]
   *
   * @param enterCommodity the enering variale
   *
   * @return exit base
   */

  Exit_base getExitBasebyPath(const ENTER_VARIABLE &enterCommodity) {
    /**
     * A lambda= beta
     *
     */
    fill(Lambda, Lambda + K + N + S, 0);

    lambda_K = Lambda;
    lambda_S = Lambda + K;
    lambda_N = Lambda + K + S;

    const vector<int> &path = enterCommodity.path;

    int commodity_primary_path_loc = primary_path_loc[enterCommodity.id];

    const vector<int> &commodity_path = paths[commodity_primary_path_loc].path;

    /**
     *  saturate link set is  empty
     *
     */
    if (0 == S) {
      /**
       * [ I_{K*K}   0       ] [ lambda_K ]  =[ beta_K ]
       * [ H         I_{N*N }] [ lambda_N ] = [ beta_N ]
       *
       *  lambda_K               = beta_K
       *  H lambda_K + lambda_N  = beta_N
       *  lambda_N               = beta_N -  H lambda_K
       *  lambda_N               = beta_N - (b_enterCommodity.id )_N
       */

      lambda_K[enterCommodity.id] = 1.0;

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        lambda_N[un_saturate_link_ids[*it]] = 1.0;
      }

      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        lambda_N[un_saturate_link_ids[*it]] -= 1.0;
      }

    } else {
      /**
       * [  I_{K*K}   B         0       ]  [ lambda_K ] = [ beta_K ]
       * [  C         D         0       ]  [ lambda_S ] = [ beta_S ]
       * [  H         F         I_{N*N} ]  [ lambda_N ] = [ beta_N ]
       *
       * lambda_S=(C beta_K-beta_S) /( CB-D )
       */
      int nrhs = 1;
      int lda = S;

      int ldb = S;
      int info = 0;

      /**
       * b=C beta_K -beta_S
       *
       */


      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        int i = saturate_link_ids[*it];

        if (i > -1) {
          lambda_S[i] = 1.0;
        }
      }

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        int i = saturate_link_ids[*it];

        if (i > -1) {
          lambda_S[i] -= 1.0;
        }
      }
      int nonzero_bata = 0;
      for (int i = 0; i < S; i++) {
        if (fabs(lambda_S[i]) > 1e-6) {
          nonzero_bata++;
        }
      }

#ifdef DEBUG

      cout << "path dimension: " << S << " nonzero elements: " << nonzero_bata
           << endl;
#endif

      /**
       * lambda_S=( C beta_K-beta_S)/( CB - D )=b/SM
       *
       */

      double t = 0;
      if (nonzero_bata > 0) {
        if (KLU == para.solver) {
          callTime(t, klusolver.solve(lambda_S));

        } else if (LAPACK == para.solver) {
          copy(SM, SM + S * S, workS);
          callTime(t, dgesv_(&S, &nrhs, workS, &lda, ipiv, lambda_S, &ldb, &info));

          if (info > 0) {
            assert(false);
            printf(
                "The diagonal element of the "
                "triangular factor of "
                "Lambda,\n");
            printf(
                "U(%i,%i) is zero, so that Lambda is "
                "singular;\n",
                info, info);
            printf(
                "the solution could not be "
                "computed.\n");
            exit(1);
          }
        } else if (SPARSE == para.solver) {
          callTime(t, sparseSolver.locSolver(lambda_S));
        } else {
          assert(false);
        }
      }
      sdata.lpsolvertime += t;


      /**
       * lambda_K=beta_K-B lambda_S
       *
       */

      lambda_K[enterCommodity.id] = 1.0;
      if (nonzero_bata > 0) {
        for (int i = 0; i < K; i++) {
          for (vector<int>::iterator it = demand_secondary_path_locs[i].begin();
               it != demand_secondary_path_locs[i].end(); it++) {
            int pindex = *it;
            int link = paths[pindex].link;
            lambda_K[i] -= lambda_S[saturate_link_ids[link]];
          }
        }
      }

      /**
       * lambda_N=beta_N- H lambda_K-F lambda_S
       *
       */

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        int i = un_saturate_link_ids[*it];

        if (i > -1) {
          lambda_N[i] = 1.0;
        }
      }

      computeLeftN(Lambda);
    }

    return computeExitVariable();
  }

  /**
   * Ax=rhs
   *
   * [  I_{K*K}   B            0       ]  [ x_K ] = [ rhs_K ]
   * [  C         D            0       ]  [ x_S ] = [ rhs_S ]
   * [  H         F            I_{N*N} ]  [ x_N ] = [ rhs_N ]
   *
   *
   * x_N=(B rhs_K- rhs_S) /( CB-D )
   *
   * [  I_{K*K}   B           0       ]  [ x_K ] = [ rhs_K ]  bandwidth
   *envery
   *demand
   * [  C         D           0       ]  [ x_S ] = [ rhs_S ]  capacity of
   *saturate
   *links
   * [  H         F           I_{N*N} ]  [ x_N ] = [ rhs_N ]  capacity of
   *other
   *links
   *rhs  =  [rhs_K  rhs_S   rhsa_N]
   *x_N  =  ( B rhs_K -rhs_S )/(CB-D)
   * @param enterCommodity the enering variale
   *
   * @return exit base
   */
  Exit_base getExitBasebyStatusLink(const ENTER_VARIABLE &exitLink) {
    /**
     *A x=rhs
     *A lambda= beta
     */

    lambda_K = Lambda;
    lambda_S = Lambda + K;
    lambda_N = Lambda + K + S;
    fill(lambda_S, lambda_S +S+ N, 0.0);
    int nrhs = 1;
    int lda = S;
    int ldb = S;
    int info = 0;

    /**
     * b= -beta_S
     *
     */

    lambda_S[saturate_link_ids[exitLink.id]] = -1.0;

#ifdef DEBUG
    cout << "link dimension: " << S << " nonzero elements: " << 1 << endl;
#endif

    /**
     * lambda_S=( C beta_K-beta_S)/( CB - D )=b/SM
     *
     */
    double t = 0;
    if (KLU == para.solver) {
      callTime(t, klusolver.solve(lambda_S));

    } else if (LAPACK == para.solver) {
      copy(SM, SM + S * S, workS);
      callTime(t, dgesv_(&S, &nrhs, workS, &lda, ipiv, lambda_S, &ldb, &info));

      if (info > 0) {
        assert(false);
        printf(
            "The diagonal element of the triangular "
            "factor of "
            "Lambda,\n");
        printf(
            "U(%i,%i) is zero, so that Lambda is "
            "singular;\n",
            info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
    } else if (SPARSE == para.solver) {
      callTime(t, sparseSolver.locSolver(lambda_S));
    }
    sdata.lpsolvertime += t;

    /**
     * lambda_K=-B lambda_S
     *
     */

    for (int i = 0; i < K; i++) {
      lambda_K[i] = 0;
      for (vector<int>::iterator it = demand_secondary_path_locs[i].begin();
           it != demand_secondary_path_locs[i].end(); it++) {
        int pindex = *it;
        int link = paths[pindex].link;
        assert(link >= 0);
        lambda_K[i] -= lambda_S[saturate_link_ids[link]];
      }
    }

    /**
     * lambda_N =  -H lambda_K-F lambda_S
     *
     */
    computeLeftN(Lambda);

    return computeExitVariable();
  }

  void pivot(ENTER_VARIABLE &entering_commodity, Exit_base &exit_base) {
    sdata.enter = entering_commodity.id;
    sdata.exit = exit_base.id;
    if (PATH_T == entering_commodity.type) {
      if (DEMAND_T == exit_base.type) {
        sdata.pivotType = NOCHANGE;
        int exit_commodity_id = exit_base.id;
        int exit_primary_pid = primary_path_loc[exit_commodity_id];

        /**
         * exit primary will been deleted from base
         * matrix
         *
         */
        deletePrimarySatuateLinks(exit_commodity_id);

        /**
         * when entering commodity and exit commodity are
         *same then replace the
         *commodity primary path with  entering path
         *
         */

        paths[exit_primary_pid].path.swap(entering_commodity.path);
        paths[exit_primary_pid].owner = entering_commodity.id;

        if (entering_commodity.id != exit_base.id) {
          /**
           * when entering commodity and exit
           *commodity are diff then replace the
           *commodity primary path with the second
           *path of the commodity which
           * crossponding saturate link cross entering
           *path and make entering path as
           *path coreesponding saturate link
           *
           */

          assert(!demand_secondary_path_locs[exit_commodity_id].empty());
          vector<int>::const_iterator it =
              demand_secondary_path_locs[exit_commodity_id].begin();

          int pid = *it;
          int link = paths[pid].link;

          assert(link >= 0);

          /**
           * from second paths to primary path
           *
           */

          demand_secondary_path_locs[exit_commodity_id].erase(it);
          primary_path_loc[exit_commodity_id] = pid;
          paths[pid].link = -1;

          demand_secondary_path_locs[entering_commodity.id].push_back(
              exit_primary_pid);

          setStatusLink(link, exit_primary_pid);
        }

        addPrimarySaturateLink(exit_commodity_id);

      } else if (STATUS_LINK == exit_base.type) {
        sdata.pivotType = NOCHANGE;
        int pid = saturate_link_path_loc[exit_base.id];
        paths[pid].path.swap(entering_commodity.path);

        /**
         * when the owner of the exit path is not owner
         * of entering path
         *
         */

        if (entering_commodity.id != paths[pid].owner) {

          deleteElem(demand_secondary_path_locs[paths[pid].owner], pid);

          demand_secondary_path_locs[entering_commodity.id].push_back(pid);
        }
        paths[pid].owner = entering_commodity.id;

      } else {
        // exit is un saturate link then from un_saturate
        // link to saturate link
        sdata.pivotType = ADDLINK;
        sdata.linkId = exit_base.id;
        if (empty_paths.empty()) {
          Path npath(entering_commodity.path, entering_commodity.id,
                     exit_base.id);

          paths.push_back(npath);

          saturate_link_ids[exit_base.id] = saturate_links.size();
          assert(find(saturate_links.begin(), saturate_links.end(),
                      exit_base.id) == saturate_links.end());

          saturate_links.push_back(exit_base.id);

          saturate_link_path_loc[exit_base.id] = paths.size() - 1;

          demand_secondary_path_locs[entering_commodity.id].push_back(
              paths.size() - 1);

        } else {
          int pid = empty_paths.back();
          empty_paths.pop_back();

          paths[pid].path.swap(entering_commodity.path);
          paths[pid].owner = entering_commodity.id;
          paths[pid].link = exit_base.id;
          assert(find(saturate_links.begin(), saturate_links.end(),
                      exit_base.id) == saturate_links.end());
          saturate_link_ids[exit_base.id] = saturate_links.size();
          saturate_links.push_back(exit_base.id);

          saturate_link_path_loc[exit_base.id] = pid;
          demand_secondary_path_locs[entering_commodity.id].push_back(pid);
        }

        int link = exit_base.id;
        addStatusLink(link);
      }

    } else {
      /**
       * entering a saturate link
       *
       */

      if (DEMAND_T == exit_base.type) {
        int enter_saturate_link = entering_commodity.id;

        int spid = saturate_link_path_loc[enter_saturate_link];
        deleteSaturateLink(enter_saturate_link);
        int exit_commodity_id = exit_base.id;

        sdata.pivotType = DELETELINK;
        sdata.linkId = exit_base.id;

        deletePrimarySatuateLinks(exit_commodity_id);
        empty_paths.push_back(primary_path_loc[exit_commodity_id]);

        if (paths[spid].owner == exit_commodity_id) {
          primary_path_loc[exit_commodity_id] = spid;

          deleteElem(demand_secondary_path_locs[exit_commodity_id], spid);
        } else {
          assert(!demand_secondary_path_locs[exit_commodity_id].empty());
          vector<int>::const_iterator it =
              demand_secondary_path_locs[exit_commodity_id].begin();

          int pid = *it;
          int link = paths[pid].link;

          assert(link >= 0);

          demand_secondary_path_locs[exit_commodity_id].erase(it);
          primary_path_loc[exit_commodity_id] = pid;

          setStatusLink(link, spid);
        }

        addPrimarySaturateLink(exit_commodity_id);

      } else if (STATUS_LINK == exit_base.type) {
        int enter_saturate_link = entering_commodity.id;

        int spid = saturate_link_path_loc[enter_saturate_link];
        deleteSaturateLink(enter_saturate_link);
        int exit_commodity_id = exit_base.id;

        sdata.pivotType = DELETELINK;
        sdata.linkId = exit_base.id;

        if (exit_base.id == entering_commodity.id) {
          empty_paths.push_back(spid);

          deleteElem(demand_secondary_path_locs[paths[spid].owner], spid);

        } else {
          int pid = saturate_link_path_loc[exit_base.id];
          empty_paths.push_back(pid);

          deleteElem(demand_secondary_path_locs[paths[pid].owner], pid);

          setStatusLink(exit_base.id, spid);
        }

      } else {
        int enter_saturate_link = entering_commodity.id;

        int spid = saturate_link_path_loc[enter_saturate_link];

        sdata.linkId = exit_base.id;

        sdata.pivotType = CHANGE_SATURATE_LINK;

        int link = exit_base.id;

        sdata.exitLink = link;

        changeSaturateLink(enter_saturate_link, link);
      }
    }
  }

  void update_edge_left_bandwith() {
    edgeLeftBandwith = update_caps;
    vector<C> allow(K, 0.0);
    /**
     * reduce primary path bandwidth
     *
     */

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      allow[paths[pid].owner] += X[i];
      const vector<int> &path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i];
      }
    }
    /**
     *
     * reduce secondary path bandwidth
     *
     */
    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      int pid = saturate_link_path_loc[link];
      allow[paths[pid].owner] += X[i + K];
      const vector<int> &path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i + K];
      }
    }

/**
 * no link overflow
 */
#ifdef DEBUG
    for (int i = 0; i < S + N; i++) {
      assert(edgeLeftBandwith[i] >= -EPS);
    }

    /**
     * every demand satifies

     */
    for (int i = 0; i < K; i++) {
      assert(allow[i] - rhs[i] >= -EPS && allow[i] - rhs[i] <= EPS);
    }

#endif
  }

  C leftBandwith(const vector<int> path) {
    if (path.empty()) {
      return 0;
    }
    C re = edgeLeftBandwith[path.front()];
    for (vector<int>::const_iterator it = path.begin(); it != path.end();
         it++) {
      re = min(re, edgeLeftBandwith[*it]);
    }

    return re;
  }

  /**
   * MU A = c
   *[ mu_K  mu_S  0 ]  [  I_{K*K}   B     0   ]  =  -[ c_K c_S 0]
   *                   [  C         D     0   ]
   *                   [  H         F  I_{N*N}]
   *
   * mu_K   + mu_S C = -c_K
   * mu_K B + mu_S D = -c_S
   *
   * mu_S( CB - D )  = c_S- c_K B
   * mu_S SM         = c_S- c_K B  = b

   * mu_K= c_K -mu_S C
   */

  void update_edge_cost() {
    int nrhs = 1;
    int lda = S;

    int ldb = S;
    int info = 0;
    /**
     * b = c_S - c_K B
     *
     */

    fill(dual_solution, dual_solution + S, 0.0);
    for (int i = 0; i < S; i++) {
      int linkid = saturate_links[i];
      int pid = saturate_link_path_loc[linkid];
      dual_solution[i] = getOrigCost(paths[pid].path);
    }

    for (int i = 0; i < S; i++) {
      int linkid = saturate_links[i];
      int pid = saturate_link_path_loc[linkid];
      int oindex = paths[pid].owner;
      int ppid = primary_path_loc[oindex];
      dual_solution[i] -= getOrigCost(paths[ppid].path);
    }

    double t = 0;
    if (KLU == para.solver) {
      callTime(t, klusolver.tsolve(dual_solution));
    } else if (LAPACK == para.solver) {
      copy(SM, SM + S * S, workS);
      callTime(t, dgesv_(&S, &nrhs, (double *)workS, &lda, ipiv, (double *)dual_solution,
                         &ldb, &info));
      if (info > 0) {
        assert(false);
        printf(
            "The diagonal element of the triangular "
            "factor of "
            "Lambda,\n");
        printf(
            "U(%i,%i) is zero, so that Lambda is "
            "singular;\n",
            info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
    } else if (SPARSE == para.solver) {
      sparseSolver.tincSolver(Mu, dual_solution);
      memcpy(Mu, dual_solution, S * sizeof(double));
    }
    sdata.lpsolvertime += t;
    update_weights = orignal_weights;


    for (int i = 0; i < S; i++) {

      update_weights[saturate_links[i]] += dual_solution[i];
    }
  }

  void printResult() {
    // verify();
    std::cout << fixed;
    if (para.info > 1) {
      sdata.using_system_time = systemTime() - sdata.start_time;
      C sobj = success_obj();
      double OBJ = computeObj();

      std::cout << systemTime() - sdata.start_time << ","
                << sobj / (TotalB + 0.0) << "," << sobj << "," << OBJ << ","
                << sdata.iterator_num << "," << sdata.empty_iterator_num << ","
                << sdata.shortestpathtime / 1000.0 << ","
                << sdata.lpsolvertime / 1000.0 << std::endl;

    } else {
      sdata.using_system_time = systemTime() - sdata.start_time;
      std::cout << "======================================" << std::endl;
      std::cout << "using time(s) :" << sdata.using_system_time << std::endl;
      std::cout << "computing shortest path use time(ms) :"
                << sdata.shortestpathtime << std::endl;
      std::cout << "solving linear equation solve use time(ms) :"
                << sdata.lpsolvertime << std::endl;
      std::cout << "iteration time: " << sdata.iterator_num << std::endl;
      std::cout << "empty iteration tiem: " << sdata.empty_iterator_num
                << std::endl;
      std::cout << "Total nonzero element: Total saturate link "
                << sdata.totalNonzero / (sdata.totalStaturateLink + 0.01)
                << std::endl;

      C sobj = success_obj();

      std::cout << "success fractional bandwidth: " << sobj << std::endl;
      std::cout << "success fractional bandwidth rat in total demand: "
                << sobj / (TotalB + 0.0) << std::endl;
    }
    if (para.info > 1) {
      std::fstream fs;
      stringstream sstr;
      if (KLU == para.solver) {
        sstr << "klusaturate" << graph.getVertex_num() << "_"
             << graph.getLink_num() << ".csv";
      } else if (LAPACK == para.solver) {
        sstr << "blassaturate" << graph.getVertex_num() << "_"
             << graph.getLink_num() << ".csv";
      } else {
        sstr << "sparsesaturate" << graph.getVertex_num() << "_"
             << graph.getLink_num() << ".csv";
      }

      fs.open(sstr.str(), std::fstream::out | std::fstream::trunc);
      for (vector<pair<int, int>>::iterator it = saturateLinkAndMatrix.begin();
           it != saturateLinkAndMatrix.end(); it++) {
        fs << it->first << "," << it->second << endl;
      }
      fs.close();
    }
  }
  bool verify(void) {
    cout << "verify" << endl;
    for (int i = 0; i < K + S; i++) {
      if (X[i] < -EPS) {
        cout << "can not have negative assigin " << endl;
      }
    }
    update_edge_left_bandwith();
    for (size_t i = 0; i < edgeLeftBandwith.size(); i++) {
      if (edgeLeftBandwith[i] < -EPS) {
        cout << "over assigin" << endl;
      }
    }
  }
};
}
}
