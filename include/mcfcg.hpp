
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

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))

#else
#include <omp.h>

#endif

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include "System.h"
#include "config.h"
#include "graphalg.hpp"
#include "klu.h"
using namespace raptor;
using namespace std;

namespace raptor {

namespace mcmcf {

struct ENTER_VARIABLE {
  ENTER_BASE_TYPE type;
  int id;
  vector<int> path;
  ENTER_VARIABLE() : type(PATH_T), id(0) {}
};

struct EXIT_VARIABLE {
  EXIT_BASE_TYPE type;
  int id;
  EXIT_VARIABLE() : type(DEMAND_T), id(0) {}
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

  int OrigLink_num;

  W Inf_weight;
  W Inf;

  vector<Demand<C>> demands;

  vector<bool> Orig_success;

  vector<W> orignal_weights;

  vector<C> orignal_caps;

  vector<W> update_weights;

  vector<C> update_caps;

  vector<C> edgeLeftBandwith;

  C TotalB;

  vector<Path> paths;  // all the paths save in this vector

  vector<int> empty_paths;  // the location which  is delete path

  vector<C> dual_solution;  // only link have dual value

  vector<double> min_commodity_cost;

  vector<int> status_links;  // contain the index of status links

  vector<int> un_status_links;

  vector<unordered_set<int>>
      demand_secondary_path_locs;  // belong to fixed demands' second paths

  unordered_map<int, unordered_set<int>>
      status_primary_path_locs;  // primary paths which
                                 // corross the status link

  vector<int> primary_path_loc;  // every demands have a primary path

  vector<int> status_link_path_loc;  //  the path corresponding status link

  Statistics_data sdata;

  vector<int> candidate_enter;

  vector<C> rhs;

  double *X;
  double *Lambda;
  double *Mu;

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

  double *mu_K;
  double *mu_S;
  double *mu_N;

  C CZERO;
  /**
   * K: number of commodity
   * S: number of saturate link
   * N: number of unsaturate link
   */

  int K, S, N; 
  double EPS;
  solverPara para;


  int thread_num;

  LU_SOLVER lu_sover;

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
      if (paths[pid].path.back() < OrigLink_num) {
        re += x_K[i];
      }
    }

    for (int i = 0; i < S; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];

      if (paths[pid].path.back() < OrigLink_num) {
        re += x_S[i];
      }
    }
    return re;
  }
  
  double computeOBJ() {
    double re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      double p_cost = 0;
      for (vector<int>::const_iterator it = paths[pid].path.begin();
           it != paths[pid].path.end(); it++) {
        p_cost += orignal_weights[*it];
      }
      re += p_cost * x_K[i];
    }

    for (int i = 0; i < S; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];

      double p_cost = 0;
      for (vector<int>::const_iterator it = paths[pid].path.begin();
           it != paths[pid].path.end(); it++) {
        p_cost += orignal_weights[*it];
      }
      re += p_cost * x_S[i];
    }
    return re;
  }
  double computeObj() {
    double re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      if (paths[pid].back() < OrigLink_num) {
        double p_cost = 0;
        for (vector<int>::const_iterator it = paths[pid].begin();
             it != paths[pid].end(); it++) {
          p_cost += orignal_weights[*it];
        }
        re += p_cost * x_K[i];
      }
    }

    for (int i = 0; i < S; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      if (paths[pid].back() < OrigLink_num) {
        double p_cost = 0;
        for (vector<int>::const_iterator it = paths[pid].begin();
             it != paths[pid].end(); it++) {
          p_cost += orignal_weights[*it];
        }
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
  /**
   * row primary
   *
   */
  void computeS() {
    if (0 == S) {
      return;
    }
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
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      const vector<int> &path = paths[pid].path;
      for(vector<int>::const_iterator lit=path.begin(); lit!=path.end(); lit++){
        vector<int>::iterator it=lower_bound(status_links.begin(), status_links.end(), *lit);
        if(it!=status_links.end()){
          int j=it-status_links.begin();
          SM[j * S + i] = -1.0;
          sdata.nzn++;
        }
      }
      // for (int j = 0; j < S; j++) {
      //   int link1 = status_links[j];
      //   if (lower_bound(path.begin(), path.end(), link1) != path.end()) {
      //     SM[j * S + i] = -1.0;
      //     sdata.nzn++;
      //   }
      // }
    }

    /**
     *  CB
     *
     */
    for (unordered_map<int, unordered_set<int>>::const_iterator it =
             status_primary_path_locs.begin();
         it != status_primary_path_locs.end(); it++) {
      int link = it->first;
      int i = getNindex(link) - K;
      for (unordered_set<int>::const_iterator pit = it->second.begin();
           pit != it->second.end(); pit++) {
        int oindex = paths[*pit].owner;
        for (unordered_set<int>::const_iterator cit =
                 demand_secondary_path_locs[oindex].begin();
             cit != demand_secondary_path_locs[oindex].end(); cit++) {
          int slink = paths[*cit].link;
          int j = getNindex(slink) - K;
          SM[i * S + j] += 1.0;
          sdata.nzn++;
        }
      }
    }
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

  void computeRHS() {
    fill(X, X + K + S + N, 0.0);
    x_K = X;
    x_S = X + K;
    x_N = X + K + S;

    int nrhs = 1;
    int lda = S;

    int ldb = S;
    int info = 0;
    /*A x = rhs
     * [  I_{K*K}     B         0     ]  [ x_K ] = [ rhs_K ]
     *bandwidth  envery
     *demand
     * [  C           D         0     ]  [ x_S ] = [ rhs_S ]
     *capacity of
     *status links
     * [  H           F       I_{N*N} ]  [ x_N ] = [ rhs_N ]
     *capacity of other
     *links
     * x_S=( C rhs_K -rhs_S )/( CB - D)
     *unordered_set SM=CB - D, b= C rhs_K -rhs_S
     *x_S = b/SM
     */
    if (S > 0) {
      fill(b, b + S, 0.0);

      for (int i = 0; i < S; i++) {
        b[i] = -rhs[K + status_links[i]];
      }

      for (int i = 0; i < S; i++) {
        int link = status_links[i];

        const unordered_set<int> &pps = status_primary_path_locs[link];
        for (unordered_set<int>::const_iterator it = pps.begin();
             it != pps.end(); it++) {
          b[i] += rhs[paths[*it].owner];
        }
      }
      copy(SM, SM + S * S, workS);

      if (KLU == lu_sover) {
        solveLP(workS, S, b);
      } else if (LAPACK == lu_sover) {
        dgesv_(&S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
               &info);
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
      } else {
        assert(false);
      }

      memcpy(x_S, b, S * sizeof(double));
    }

    /**
     * x_K=rhs_K-B x_S
     *
     */
    for (int i = 0; i < K; i++) {
      x_K[i] = rhs[i];
    }

    for (int i = 0; i < K; i++) {
      const unordered_set<int> &pathindices = demand_secondary_path_locs[i];
      for (unordered_set<int>::const_iterator it = pathindices.begin();
           it != pathindices.end(); it++) {
        x_K[i] -= x_S[getNindex(paths[*it].link) - K];
      }
    }
  }

  /**
   *
   * @brief choose candiate leaving base
   *
   * @return
   */
  EXIT_VARIABLE computeLeavingVariable() {
    /**
     *  choose enter base as i=argmin_{ lambda[ i]>0} x[ i ]/lambda[ i ]
     *
     */

    int reK = -1;
    double minK = numeric_limits<double>::max();
    double temp;
    int i = 0;
    while (i < K) {
      if (lambda_K[i] > EPS) {
        reK = i;
        minK = x_K[i] / lambda_K[i];
        break;
      }
      i++;
    }
    while (i < K) {
      if (lambda_K[i] > EPS) {
        temp = x_K[i] / lambda_K[i];
        if (temp < minK) {
          reK = i;
          minK = temp;
        }
      }
      i++;
    }

    int reN = -1;
    double minN = numeric_limits<double>::max();
    i = 0;
    while (i < S) {
      if (lambda_S[i] > EPS) {
        reN = i;
        minN = x_S[i] / lambda_S[i];
        break;
      }
      i++;
    }

    while (i < S) {
      if (lambda_S[i] > EPS) {
        temp = x_S[i] / lambda_S[i];
        if (temp < minN) {
          reN = i;
          minN = temp;
        }
      }
      i++;
    }

    int reJ = -1;
    double minJ = numeric_limits<double>::max();

    i = 0;
    while (i < N) {
      if (lambda_N[i] > EPS) {
        reJ = i;
        minJ = x_N[i] / lambda_N[i];
        break;
      }
      i++;
    }

    while (i < N) {
      if (lambda_N[i] > EPS) {
        temp = x_N[i] / lambda_N[i];
        if (temp < minJ) {
          reJ = i;
          minJ = temp;
        }
      }
      i++;
    }

    EXIT_VARIABLE re;

    if (minK <= minN && minK <= minJ) {
      re.id = reK;
      re.type = DEMAND_T;
      if (minK <= EPS) sdata.empty_iterator_num++;
      return re;
    }
    if (minN <= minJ) {
      re.id = status_links[reN];
      re.type = STATUS_LINK;
      if (minN <= EPS) {
        sdata.empty_iterator_num++;
      }
      return re;
    }
    re.id = un_status_links[reJ];
    re.type = OTHER_LINK;
    if (minJ <= EPS) {
      sdata.empty_iterator_num++;
    }
    return re;
  }

 public:
  CG(const G &g, const vector<W> &ws, const vector<C> &caps,
     const vector<Demand<C>> &ds)
      : demands(ds),  thread_num(1), lu_sover(KLU) {
    orignal_weights.resize(ws.size());
    orignal_caps.resize(ws.size());
    vector<int> srcs, snks;
    linkMap.resize(ws.size());
    /**
     *  Reorder links' id to improve
     * graph get adjacent nodes speed
     *
     */

    for (size_t v = 0; v < g.getVertex_num(); v++) {
      int degree = g.getOutDegree(v);
      for (int i = 0; i < degree; i++) {
        int link = g.getAdj(v, i);
        linkMap[srcs.size()] = link;
        orignal_caps[srcs.size()]=caps[link];
        orignal_weights[srcs.size()] = ws[link];
        int snk;
        g.findRhs(link, v, snk);
        srcs.push_back(v);
        snks.push_back(snk);
      }
    }
    graph.initial(srcs, snks);

    double max_w = 0;
    for (int i = 0; i < OrigLink_num; i++) {
      if (max_w < orignal_weights[i]) {
        max_w = orignal_weights[i];
      }
    }

    Inf_weight = max_w * (graph.getVertex_num() - 1) + 1;
    Inf = Inf_weight;

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < 1000; i++) {
      thread_num = omp_get_num_threads();
    }

#endif  // _OPENMP


    OrigLink_num = graph.getLink_num();

    status_link_path_loc.resize(OrigLink_num+1);

    CZERO = ((C)1e-6);
    EPS = ((double)1e-4);

    K = demands.size();
    demand_secondary_path_locs.resize(K);

    S = 0;
    N = OrigLink_num + 1;

    b = new double[N];
    X = new double[K + N];
    Lambda = new double[K + N];
    Mu=new double[K + N];

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

    if (minNonzero > 1) {
      for (size_t i = 0; i < update_caps.size(); i++) {
        update_caps[i] += 0.0001 * (rand() % 1000) / 1000.0;
      }

    } else {
      for (size_t i = 0; i < update_caps.size(); i++) {
        update_caps[i] += minNonzero * (rand() % 1000) / 1000.0;
      }
    }

    orignal_weights.push_back(Inf_weight);
    update_caps.push_back(TotalB);
  }
  ~CG() {
    if (NULL != Lambda) {
      delete[] Lambda;
      Lambda = NULL;
    }

    if (NULL != X) {
      delete[] X;
      X = NULL;
    }
    if(NULL!= Mu){
      delete [] Mu;
      Mu=NULL;
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
  void setInfo(const int level) { para.info = level; }

  void setLUSOLVER(LU_SOLVER s) { lu_sover = s; }

  bool is_status_link(const int link) const {
    return lower_bound(status_links.begin(), status_links.end(), link) !=
           status_links.end();
  }

  void setStatusLink(const int link, const int pid) {
    paths[pid].link = link;
    status_link_path_loc[link] = pid;
  }

  void deletePrimaryPath(const int commodityID) {
    int primary_pid = primary_path_loc[commodityID];
    const vector<int> &orig_path = paths[primary_pid].path;
    for (vector<int>::const_iterator it = orig_path.begin();
         it != orig_path.end(); it++) {
      if (lower_bound(status_links.begin(), status_links.end(), *it) !=
          status_links.end()) {
        status_primary_path_locs[*it].erase(primary_pid);
      }
    }
  }

  void addPrimaryPath(const int commodityID) {
    int pid = primary_path_loc[commodityID];
    const vector<int> &new_path = paths[pid].path;
    for (vector<int>::const_iterator it = new_path.begin();
         it != new_path.end(); it++) {
      if (lower_bound(status_links.begin(), status_links.end(), *it) !=
          status_links.end()) {
        status_primary_path_locs[*it].insert(pid);
      }
    }
  }

  void addStatusLink(const int link) {
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      if (lower_bound(paths[pid].path.begin(), paths[pid].path.end(), link) !=
          paths[pid].path.end()) {
        status_primary_path_locs[link].insert(pid);
      }
    }
  }
  void deleteStatusLink(const int link) {
    vector<int>::iterator it =
        lower_bound(status_links.begin(), status_links.end(), link);
    assert(it != status_links.end());

    status_links.erase(it);

    int spid = status_link_path_loc[link];
    paths[spid].link = -1;

    status_primary_path_locs.erase(link);
  }
  /**
   *
   *
   * @param i the id of un_statrus link
   *
   * @return  the index of the un_status link in simplex matrix
   */
  int getJIndex(int i) const {
    assert(lower_bound(un_status_links.begin(), un_status_links.end(), i) !=
           un_status_links.end());
    vector<int>::const_iterator it =
        lower_bound(un_status_links.begin(), un_status_links.end(), i);
    return (it - un_status_links.begin()) + K + S;
  }
  /**
   *
   *
   * @param i the id of status link
   *
   * @return the index of the status link in smplex matrix
   */
  int getNindex(int i) const {
    assert(lower_bound(status_links.begin(), status_links.end(), i) !=
           status_links.end());
    vector<int>::const_iterator it =
        lower_bound(status_links.begin(), status_links.end(), i);
    return (it - status_links.begin()) + K;
  }

  double getOrigCost(const vector<int> &path) const {
    double re = 0.0;
    for (vector<int>::const_iterator it = path.begin(); it != path.end();
         it++) {
      re += orignal_weights[*it];
    }
    return re;
  }

  void computIndexofLinks() {
    stable_sort(status_links.begin(), status_links.end());
    un_status_links.clear();

    size_t i = 0;
    for (int j = 0; j < S + N; j++) {
      if (i >= status_links.size()) {
        un_status_links.push_back(j);
      } else {
        if (j < status_links[i]) {
          un_status_links.push_back(j);
        } else {
          i++;
        }
      }
    }
  }
  /**
   *
   *
   * @param M  column master
   * @param n
   * @param b
   *
   * @return
   */
  bool solveLP(const double *M, const int n, double *b) {
    vector<int> AAp;
    vector<int> AAi;
    vector<double> AAx;
    int nz = 0;
    AAp.push_back(0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (fabs(M[i * n + j]) > EPS) {
          nz++;
          AAi.push_back(j);
          AAx.push_back(M[i * n + j]);
        }
      }
      AAp.push_back(nz);
    }

    klu_symbolic *Symbolic;
    klu_numeric *Numeric;
    klu_common Common;

    int *Ap = &AAp[0];
    int *Ai = &AAi[0];
    double *Ax = &AAx[0];
    klu_defaults(&Common);
    Symbolic = klu_analyze(n, Ap, Ai, &Common);
    Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);
    klu_solve(Symbolic, Numeric, n, 1, b, &Common);
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);

    return true;
  }

  bool solve() {
    sdata.start_time = systemTime();
    initial_solution();

    iteration();

    if (para.info > 0) {
      printResult();
    }
    return true;
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
      for (int i = 0; i < OrigLink_num; i++) {
        if (orignal_caps[i] < bw) {
          ws[i] = 2 * Inf_weight;
        }
      }

      inc_ksp::yen_ksp<G, vector<W>, W> ksp_g(graph, orignal_weights,
                                              Inf_weight);
      inc_ksp::yen_next_path<G, vector<W>, W> ksp = ksp_g.next_path(src, snk);
      vector<int> path;
      for (int j = 0; j < k && ksp.next_path(path); j++) {
        paths[i].push_back(path);
      }
    }
  }

  /**
   * Set a initial routing path for every demand, the only constraint the
   * initial solution must satisfies is that link capacity.
   *
   */

  void initial_solution() {
    paths.resize(K);

    primary_path_loc.resize(K, 0);

    Orig_success.resize(K, true);
    vector<W> ws(orignal_weights.size(), 1);

/**
 *  check whether there is path from  demand source to demand target
 *
 */

#pragma omp parallel for
 
    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      vector<int> path;
      if (!bidijkstra_shortest_path(graph, ws, src, snk, path,
                                    graph.getLink_num() + 1)) {
        Orig_success[i] = false;
      }
    }

    vector<C> temp_cap(orignal_caps);
    vector<bool> success(K, false);

    /**
     * Direcltly setup demand one by one
     *
     */

    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      vector<W> ws = orignal_weights;

      for (int j = 0; j < OrigLink_num; j++) {
        if (temp_cap[j] < bw) {
          ws[j] = Inf_weight;
        }
      }

      vector<int> path;
      if (bidijkstra_shortest_path(graph, ws, src, snk, path, Inf_weight)) {
        success[i] = true;

        for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
          temp_cap[*it] -= bw;
        }
        sort(path.begin(), path.end());
        paths[i].path = path;
        paths[i].owner = i;
        primary_path_loc[i] = i;
      }
    }

    for (int i = 0; i < K; i++) {
      if (!success[i]) {
        vector<int> path;
        path.push_back(OrigLink_num);
        paths[i].path = path;
        paths[i].owner = i;
        primary_path_loc[i] = i;
        success[i] = true;
      }
    }



    update_weights = orignal_weights;

    min_commodity_cost.resize(K, 2 * Inf_weight);
#pragma omp parallel for
    for (int i = 0; i < K; i++) {
      vector<int> path;
      int src = demands[i].src;
      int snk = demands[i].snk;
      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                                   2 * Inf_weight)) {
        min_commodity_cost[i] = path_cost(update_weights, path, (W)0.0);
      }
    }


    for (int i = 0; i < N; i++) {
      un_status_links.push_back(i);
    }

    rhs.resize(K + OrigLink_num+1, (C)0.0);
    for (int i = 0; i < K; i++) {
      rhs[i] = demands[i].bandwidth;
    }

    for (size_t i = 0; i < update_caps.size(); i++) {
      rhs[i + K] = update_caps[i];
    }

    dual_solution.resize(N, 0);


  }

  void iteration() {
    while (true) {
      sdata.iterator_num++;
      if (para.info > 0) {
        if (sdata.iterator_num % para.perIterationPrint == 0) {
          sdata.using_system_time = systemTime() - sdata.start_time;
          C sobj = success_obj();
          std::cout << "using time(s) :" << sdata.using_system_time
                    << std::endl;
          std::cout << sdata.iterator_num << " status link num: " << S
                    << " objvalue: " << computeOBJ()
                    << " success fractional bw: " << sobj
                    << " success rat: " << sobj / (TotalB + 0.01)
                    << " the obj gap from opt is: " << sdata.estimee_opt_diff
                    << std::endl;
          std::cout << "nonzero matrix values: " << sdata.nzn
                    << "   SLU non-zeros " << sdata.snzn << std::endl;
          std::cout << sdata.etype << " enter: " << sdata.enter << "    "
                    << sdata.exitt << " exit: " << sdata.exit << std::endl;
        }
      }
      if (sdata.iterator_num > para.maxIterationNum) {
        return;
      }

      /**
       *  entering variable choose
       *
       */
      update_edge_left_bandwith();

      ENTER_VARIABLE entering_commodity = chooseEnteringVariable();
      if (entering_commodity.id < 0) {
        return;
      }

      EXIT_VARIABLE leaving_base;

      if (PATH_T == entering_commodity.type) {
        /**
         *  leaving base  choose
         *
         */
        leaving_base = getLeavingBasebyPath(entering_commodity);
      } else {
        leaving_base = getLeavingBasebyStatusLink(entering_commodity);
      }

      pivot(entering_commodity, leaving_base);

      S = status_links.size();
      N = OrigLink_num + 1 - S;
      computIndexofLinks();

      computeS();

      update_edge_cost();

      /**
       * column primary
       *
       */

      transposeS();

      computeRHS();
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
  ENTER_VARIABLE chooseEnteringVariable() {
    ENTER_VARIABLE enter_variable;

    enter_variable.type = PATH_T;
    enter_variable.id = -1;
    /**
     *  check whether saturate link dual value is negative
     *
     */
    double max_diff = -EPS;
    for (int i = 0; i < S; i++) {
      int link = status_links[i];
      if (dual_solution[link] < max_diff) {
        max_diff = dual_solution[link];
        enter_variable.id = link;
        enter_variable.type = LINK_T;
      }
    }
    
    /**
     * If there is a dual value of  saturate link is negative then this link is a entering variable
     */

    if (enter_variable.id >= 0) {
      return enter_variable;
    }

    /**
     * fast way to choose entering variable as a path
     */
    for (size_t i = 0; i < candidate_enter.size(); i++){

      int id = candidate_enter[i];
      int src = demands[id].src;
      int snk = demands[id].snk;
      vector<int> path;

      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[id]].path, (W)0.0);

      bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                               Inf_weight);
      W new_cost = path_cost(update_weights, path, (W)0.0);
      W diff=old_cost - new_cost;

      if (diff > 10000 * EPS &&  diff> sdata.objSpeed ) {
        candidate_enter.erase(candidate_enter.begin(), candidate_enter.begin()+i+1 );
        enter_variable.id = id;
        enter_variable.type = PATH_T;
        enter_variable.path = path;
        sort(enter_variable.path.begin(), enter_variable.path.end());
        sdata.objSpeed =(para.objSpeedUpdateRat)*sdata.objSpeed +(1-(para.objSpeedUpdateRat))*diff;
        return enter_variable;
      }
    }
    candidate_enter.clear();

    sdata.estimee_opt_diff = 0;
    vector<double> opt_gap(thread_num, 0);
    vector<double> max_gaps(thread_num, EPS);
    vector<double> max_diffs(thread_num, EPS);
    vector<vector<int>> candidate(thread_num);

    vector<vector<int>> path(thread_num);
    vector<ENTER_VARIABLE> enter_variables(thread_num);
    
    for (int i = 0; i < thread_num; i++) {
      enter_variables[i].type = PATH_T;
      enter_variables[i].id = -1;
    }
#pragma omp parallel for
    for (int i = 0; i < K; i++) {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;

#endif  // (_OPENMP)
      if (!Orig_success[i]) {
        continue;
      }
      int src = demands[i].src;
      int snk = demands[i].snk;
      max_diff = EPS;
      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]].path, (W)0.0);
      /**
       * If the possible biggest reduce objective value is less than exist found one then continue
       * 
       */

      if ((old_cost - min_commodity_cost[i]) * demands[i].bandwidth <
          max_gaps[tid]) {
        continue;
      }

      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path[tid],
                                   Inf_weight)) {
        W new_cost = path_cost(update_weights, path[tid], (W)0.0);

        if (new_cost > Inf) {
          path[tid].clear();
          path[tid].push_back(OrigLink_num);
          new_cost = Inf;
        }

        W temp_diff = (old_cost - new_cost);

        if (temp_diff > EPS) {
          if (temp_diff > max_diffs[tid]) {
            max_diffs[tid] = temp_diff;
          }
          opt_gap[tid] += temp_diff * demands[i].bandwidth;

          temp_diff *= leftBandwith(path[tid]) + EPS;

          if (temp_diff > max_gaps[tid]) {
            if (temp_diff > 10000 * EPS) {
              candidate[tid].push_back(i);
            }

            max_gaps[tid] = temp_diff;

            enter_variables[tid].id = i;
            enter_variables[tid].path = path[tid];
          }
        }
      }
    }

    int chid = 0;
    double max_gap = max_gaps[0];
    candidate_enter = candidate[0];
    for (int i = 1; i < thread_num; i++) {
      candidate_enter.insert(candidate_enter.end(), candidate[i].begin(),
                             candidate[i].end());
      if (max_gaps[i] > max_gap) {
        max_gap = max_gaps[i];
        chid = i;
      }
    }

    for (int i = 0; i < thread_num; i++) {
      sdata.estimee_opt_diff += opt_gap[i];
    }

    sort(enter_variables[chid].path.begin(), enter_variables[chid].path.end());

    sdata.objSpeed =(para.objSpeedUpdateRat)*sdata.objSpeed +(1-(para.objSpeedUpdateRat))*max_diffs[chid];
    return enter_variables[chid];
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
   * @return leaving base
   */

  EXIT_VARIABLE getLeavingBasebyPath(const ENTER_VARIABLE &enterCommodity) {
    /**
     * A X=rhs
     *
     */

    fill(X, X + K + S + N, 0.0);
    x_K = X;
    x_S = X + K;
    x_N = X + K + S;

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
      fill(Lambda, Lambda+K+N, 0);
      
      lambda_K[enterCommodity.id] = 1.0;

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        lambda_N[getJIndex(*it) - S - K] = 1.0;
      }

      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        lambda_N[getJIndex(*it) - S - K] -= 1.0;
      }


      x_K = X;
      x_N = X + K;
      /**
       * [ I_{K*K}   0       ] [ x_K ]  =[ rhs_K ]
       * [ H         I_{N*N }] [ x_N ] = [ rhs_N ]
       *
       *  x_K = rhs_K
       *  H rhs_K + x_N = rhs_N
       *  x_N = rhs_N -H x_K
       */

      for (int i = 0; i < K; i++) {
        x_K[i] = rhs[i];
      }

      /**
       * x_N = rhs_N -H x_K
       *
       */

      for (int i = 0; i < N; i++) {
        x_N[i] = edgeLeftBandwith[un_status_links[i]];
      }

    } else {
      /**
       * [  I_{K*K}   B         0       ]  [ lambda_K ] = [ beta_K ]
       * [  C         D         0       ]  [ lambda_S ] = [ beta_S ]
       * [  H         F         I_{N*N} ]  [ lambda_N ] = [ beta_N ]
       *
       * lambda_S=(B beta_K-beta_S) /( CB-D )
       */
      int nrhs = 1;
      int lda = S;

      int ldb = S;
      int info = 0;

      /**
       * b=B beta_K -beta_S
       *
       */

      fill(b, b + S, 0.0);
      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        vector<int>::iterator fid =
            lower_bound(status_links.begin(), status_links.end(), *it);
        if (fid != status_links.end()) {
          b[fid - status_links.begin()] = 1.0;
        }
      }

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        vector<int>::iterator fid =
            lower_bound(status_links.begin(), status_links.end(), *it);
        if (fid != status_links.end()) {
          b[fid - status_links.begin()] -= 1.0;
        }
      }
      /**
       * lambda_S=( C beta_K-beta_S)/( CB - D )=b/SM
       *
       */

      if (KLU == lu_sover) {
        copy(SM, SM + S * S, workS);
        solveLP(workS, S, b);
      } else if (LAPACK == lu_sover) {
        char c = 'S';
        dgetrs_(&c, &S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
                &info);

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
      } else {
        assert(false);
      }

      memcpy(lambda_S, b, S * sizeof(double));

      /**
       * lambda_K=beta_K-Lambda lambda_S
       *
       */

      lambda_K[enterCommodity.id] = 1.0;
      for (int i = 0; i < K; i++) {
        for (unordered_set<int>::iterator it =
                 demand_secondary_path_locs[i].begin();
             it != demand_secondary_path_locs[i].end(); it++) {
          int pindex = *it;
          int link = paths[pindex].link;
          lambda_K[i] -= lambda_S[getNindex(link) - K];
        }
      }
      
      /**
       * lambda_N=beta_N- H lambda_K-F lambda_S
       *
       */

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        if (!is_status_link(*it)) lambda_N[getJIndex(*it) - S - K] = 1.0;
      }

      vector<int> left_AK;
      for (int i = 0; i < K; i++) {
        if (fabs(lambda_K[i]) > EPS) {
          left_AK.push_back(i);
        }
      }

      for (size_t i = 0; i < un_status_links.size(); i++) {
        for (size_t k = 0; k < left_AK.size(); k++) {
          int pid = left_AK[k];
          int ppid = primary_path_loc[pid];
          if (lower_bound(paths[ppid].path.begin(), paths[ppid].path.end(),
                          un_status_links[i]) != paths[ppid].path.end())
            lambda_N[i] -= lambda_K[pid];
        }
      }

      vector<int> left_AN;

      for (int i = 0; i < S; i++) {
        if (fabs(lambda_S[i]) > EPS) {
          left_AN.push_back(i);
        }
      }

      for (size_t i = 0; i < un_status_links.size(); i++) {
        for (size_t k = 0; k < left_AN.size(); k++) {
          int pid = left_AN[k];
          int ppid = status_link_path_loc[status_links[pid]];
          if (lower_bound(paths[ppid].path.begin(), paths[ppid].path.end(),
                          un_status_links[i]) != paths[ppid].path.end())
            lambda_N[i] -= lambda_S[pid];
        }
      }
    }

    return computeLeavingVariable();
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
   * @return leaving base
   */
  EXIT_VARIABLE getLeavingBasebyStatusLink(const ENTER_VARIABLE &enterLink) {
    /**
     *A x=rhs
     *A lambda= beta
     */

    fill(Lambda, Lambda + K + S + N, 0.0);
    lambda_K = Lambda;
    lambda_S = Lambda + K;
    lambda_N = Lambda + K + S;

    int nrhs = 1;
    int lda = S;

    int ldb = S;
    int info = 0;

    /**
     * b= -beta_S
     *
     */
    fill(b, b + S, 0.0);

    b[getNindex(enterLink.id) - K] = -1.0;


    /**
     * lambda_S=( C beta_K-beta_S)/( CB - D )=b/SM
     *
     */
    if (KLU == lu_sover) {
      copy(SM, SM + S * S, workS);
      solveLP(workS, S, b);
    } else if (LAPACK == lu_sover) {
      char c = 'S';
      dgetrs_(&c, &S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
              &info);
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
    }

    memcpy(lambda_S, b, S * sizeof(double));

    /**
     * lambda_K=-B lambda_S
     *
     */

    for (int i = 0; i < K; i++) {
      for (unordered_set<int>::iterator it = demand_secondary_path_locs[i].begin();
           it != demand_secondary_path_locs[i].end(); it++) {
        int pindex = *it;
        int link = paths[pindex].link;
        assert(link >= 0);
        lambda_K[i] -= lambda_S[getNindex(link) - K];
      }
    }

    /**
     * lambda_N =  -H lambda_K-F lambda_S
     *
     */

    vector<int> left_AK;
    for (int i = 0; i < K; i++) {
      if (fabs(lambda_K[i]) > EPS) {
        left_AK.push_back(i);
      }
    }

    for (size_t i = 0; i < un_status_links.size(); i++) {
      for (size_t k = 0; k < left_AK.size(); k++) {
        int pid = left_AK[k];
        int ppid = primary_path_loc[pid];

        if (!paths[ppid].path.empty() &&
            lower_bound(paths[ppid].path.begin(), paths[ppid].path.end(),
                        un_status_links[i]) != paths[ppid].path.end())
          lambda_N[i] -= lambda_K[pid];
      }
    }

    vector<int> left_AN;

    for (int i = 0; i < S; i++) {
      if (fabs(lambda_S[i]) > EPS) {
        left_AN.push_back(i);
      }
    }

    for (size_t i = 0; i < un_status_links.size(); i++) {
      for (size_t k = 0; k < left_AN.size(); k++) {
        int pid = left_AN[k];
        int ppid = status_link_path_loc[status_links[pid]];
        if (!paths[ppid].path.empty() &&
            lower_bound(paths[ppid].path.begin(), paths[ppid].path.end(),
                        un_status_links[i]) != paths[ppid].path.end())
          lambda_N[i] -= lambda_S[pid];
      }
    }

    return computeLeavingVariable();
  }

  void pivot(ENTER_VARIABLE &entering_commodity, EXIT_VARIABLE &leaving_base) {
    if (PATH_T == entering_commodity.type) {
      if (DEMAND_T == leaving_base.type) {
        int exit_commodity_id = leaving_base.id;
        int exit_primary_pid = primary_path_loc[exit_commodity_id];

        /**
         * exit primary will been deleted from base
         * matrix
         *
         */
        deletePrimaryPath(exit_commodity_id);

        /**
         * when enter commodity and exit commodity are
         *same then replace the
         *commodity primary path with  enter path
         *
         */

        paths[exit_primary_pid].path = entering_commodity.path;
        paths[exit_primary_pid].owner = entering_commodity.id;

        if (entering_commodity.id != leaving_base.id) {
          /**
           * when enter commodity and exit
           *commodity are diff then replace the
           *commodity primary path with the second
           *path of the commodity which
           * crossponding status link cross enter
           *path and make enter path as
           *path coreesponding status link
           *
           */

          assert(!demand_secondary_path_locs[exit_commodity_id].empty());
          unordered_set<int>::const_iterator it =
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

          demand_secondary_path_locs[entering_commodity.id].insert(
              exit_primary_pid);

          setStatusLink(link, exit_primary_pid);
        }

        addPrimaryPath(exit_commodity_id);

      } else if (STATUS_LINK == leaving_base.type) {
        int pid = status_link_path_loc[leaving_base.id];
        paths[pid].path = entering_commodity.path;

        /**
         * when the owner of the exit path is not owner
         * of enter path
         *
         */

        if (entering_commodity.id != paths[pid].owner) {
          demand_secondary_path_locs[paths[pid].owner].erase(pid);
          demand_secondary_path_locs[entering_commodity.id].insert(pid);
        }
        paths[pid].owner = entering_commodity.id;

      } else {
        // exit is un status link then from un_status
        // link to status link
        if (empty_paths.empty()) {
          Path npath(entering_commodity.path, entering_commodity.id,
                     leaving_base.id);

          paths.push_back(npath);

          status_links.push_back(leaving_base.id);
          status_link_path_loc[leaving_base.id] = paths.size() - 1;

          demand_secondary_path_locs[entering_commodity.id].insert(paths.size() -
                                                                1);

        } else {
          int pid = empty_paths.back();
          empty_paths.pop_back();

          paths[pid].path = entering_commodity.path;
          paths[pid].owner = entering_commodity.id;
          paths[pid].link = leaving_base.id;

          status_links.push_back(leaving_base.id);

          status_link_path_loc[leaving_base.id] = pid;
          demand_secondary_path_locs[entering_commodity.id].insert(pid);
        }

        int link = leaving_base.id;
        addStatusLink(link);
      }

    } else {
      /**
       * enter a status link
       *
       */
      int enter_status_link = entering_commodity.id;

      int spid = status_link_path_loc[enter_status_link];
      deleteStatusLink(enter_status_link);

      if (DEMAND_T == leaving_base.type) {
        int exit_commodity_id = leaving_base.id;

        deletePrimaryPath(exit_commodity_id);

        if (paths[spid].owner == leaving_base.id) {
          primary_path_loc[leaving_base.id] = spid;
          demand_secondary_path_locs[leaving_base.id].erase(spid);
        } else {
          assert(!demand_secondary_path_locs[leaving_base.id].empty());
          unordered_set<int>::const_iterator it =
              demand_secondary_path_locs[leaving_base.id].begin();

          int pid = *it;
          int link = paths[pid].link;

          assert(link >= 0);

          demand_secondary_path_locs[leaving_base.id].erase(it);
          primary_path_loc[leaving_base.id] = pid;

          setStatusLink(link, spid);
        }

        addPrimaryPath(exit_commodity_id);
      } else if (STATUS_LINK == leaving_base.type) {
        if (leaving_base.id == entering_commodity.id) {
          empty_paths.push_back(spid);
          demand_secondary_path_locs[paths[spid].owner].erase(spid);
          paths[spid].owner = -1;
        } else {
          int pid = status_link_path_loc[leaving_base.id];
          empty_paths.push_back(pid);
          demand_secondary_path_locs[paths[pid].owner].erase(pid);

          paths[pid].owner = -1;
          paths[pid].link = -1;

          setStatusLink(leaving_base.id, spid);
        }

      } else {
        status_links.push_back(leaving_base.id);

        int link = leaving_base.id;

        setStatusLink(link, spid);

        addStatusLink(link);
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
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      allow[paths[pid].owner] += X[i + K];
      const vector<int> &path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i + K];
      }
    }
    /**
     * no link over is  used
     */
#ifdef DEBUG
    for (int i = 0; i < S + N; i++) {
      assert(edgeLeftBandwith[i] >= -EPS);
    }
    
    /**
     * every demand satifies

     */
    for (int i = 0; i < K; i++) {
      assert(allow[i] - demands[i].bandwidth >= -EPS &&
             allow[i] - demands[i].bandwidth <= EPS);
    }

    x_N = X  + K+ S;
    for (int i = 0; i < N; i++) {
      int link = un_status_links[i];
      x_N[i] = edgeLeftBandwith[link];
      
      assert(x_N[i] >= -EPS);
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
   *[ mu_K  mu_S  0 ]  [  I_{K*K}   B     0   ]  =  [ c_K c_N 0]
   *                   [  C         D     0   ]
   *                   [  H         F  I_{N*N}]
   *
   * mu_K   + mu_S C = c_K
   * mu_K B + mu_S D = c_N
   *
   * mu_S( CB - D ) = c_K B - c_N
   * mu_S SM = c_K B -c_N = b
   
   * mu_K= c_K -mu_S C
   */
  
  void update_edge_cost() {
    
    int nrhs = 1;
    int lda = S;

    int ldb = S;
    int info = 0;
    /**
     * b = c_K B -c_N 
     * 
     */

    fill(b, b + S, 0.0);
    for (int i = 0; i < S; i++) {
      int linkid = status_links[i];
      int pid = status_link_path_loc[linkid];
      b[i] = -getOrigCost(paths[pid].path);
    }

    fill(Mu, Mu + K, 0.0);
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      Mu[i] = getOrigCost(paths[pid].path);
    }
    for (int i = 0; i < S; i++) {
      int linkid = status_links[i];
      int pid = status_link_path_loc[linkid];
      int oindex = paths[pid].owner;
      b[i] += Mu[oindex];
    }

    copy(SM, SM + S * S, workS);
    if (KLU == lu_sover) {
      solveLP(workS, S, b);
    } else if (LAPACK == lu_sover) {
      dgesv_(&S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb, &info);
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
    }

    update_weights = orignal_weights;

    fill(dual_solution.begin(), dual_solution.end(), 0.0);

    for (int i = 0; i < S; i++) {
      dual_solution[status_links[i]] = b[i];
      update_weights[status_links[i]] += b[i];
    }
  }

  void printResult() {
    sdata.using_system_time = systemTime() - sdata.start_time;
    std::cout << "using time(s) :" << sdata.using_system_time << std::endl;
    std::cout << "iteration time: " << sdata.iterator_num << std::endl;
    std::cout << "empty iteration tiem: " << sdata.empty_iterator_num
              << std::endl;

    C sobj = success_obj();

    std::cout << "success fractional bandwidth: " << sobj << std::endl;
    std::cout << "success fractional bandwidth rat in total demand: "
              << sobj / (TotalB + 0.0) << std::endl;
  }
};
}
}
