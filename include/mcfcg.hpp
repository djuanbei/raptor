
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

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include "System.h"
#include "graphalg.hpp"
#include "config.h"
#include "klu.h"
using namespace raptor;
using namespace std;


namespace raptor{

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

  int origLink_num;

  W inf_weight;
  W Inf;

  vector<Demand<C> > demands;

  vector<bool> orig_success;

  vector<W> orignal_weights;
  


  vector<C> orignal_caps;

  vector<W> update_weights;

  vector<C> update_caps;

  vector<C> edgeLeftBandwith;

  C totalB;

  vector<Path> paths;  // all the paths save in this vector

  vector<int> empty_paths;  // the location which  is delete path

  vector<C> dual_solution;  // only link have dual value

  vector<double> min_commodity_cost;

  vector<int> status_links;  // contain the index of status links

  vector<int> un_status_links;

  vector<unordered_set<int>>
      demand_second_path_locs;  // belong to fixed demands' second paths

  unordered_map<int, unordered_set<int>>
      status_primary_path_locs;  // primary paths which
                                 // corross the status link

  vector<int> primary_path_loc;  // every demands have a primary path

  vector<int> status_link_path_loc;  //  the path corresponding status link

  Statistics_data sdata;

  vector<int> candidate_enter;

  vector<C> rhs;
  double *A;
  double *X;

  int *ipiv;
  double *b;

  double *S;
  double *workS;
  int S_maxdim;

  double *a_K;
  double *a_N;
  double *a_J;

  double *x_K;
  double *x_N;
  double *x_J;

  double *y_K;
  double *y_N;
  double *y_J;

  C CZERO;
  int K, N, J;
  double EPS;

  int info;

  int thread_num;

  LU_SOLVER lu_sover;

  void allocateS(const int N) {
    if (S_maxdim >= N){
      return;
    }
    if (NULL != S) {
      delete[] S;
      S = NULL;
      delete[] workS;
    }
    S_maxdim = 1.3 * N + 1;
    S = new double[S_maxdim * S_maxdim];
    workS = new double[S_maxdim * S_maxdim];
  }
  C success_obj() {
    C re = 0;
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      if (paths[pid].path.back() < origLink_num){
        re += x_K[i];
      }
    }

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];

      if (paths[pid].path.back() < origLink_num){
        re += x_N[i];
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

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];

      double p_cost = 0;
      for (vector<int>::const_iterator it = paths[pid].path.begin();
           it != paths[pid].path.end(); it++) {
        p_cost += orignal_weights[*it];
      }
      re += p_cost * x_N[i];
    }
    return re;
  }
  double computeObj() {
    double re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      if (paths[pid].back() < origLink_num) {
        double p_cost = 0;
        for (vector<int>::const_iterator it = paths[pid].begin();
             it != paths[pid].end(); it++) {
          p_cost += orignal_weights[*it];
        }
        re += p_cost * x_K[i];
      }
    }

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      if (paths[pid].back() < origLink_num) {
        double p_cost = 0;
        for (vector<int>::const_iterator it = paths[pid].begin();
             it != paths[pid].end(); it++) {
          p_cost += orignal_weights[*it];
        }
        re += p_cost * x_N[i];
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
    if (0 == N){
      return;
    }
    /**
     * S=BA-E
     *
     */
    sdata.nzn = 0;
    if (N > S_maxdim){
      allocateS(N);
    }
    fill(S, S + N * N, 0.0);

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      const vector<int> &path = paths[pid].path;
      for (int j = 0; j < N; j++) {
        int link1 = status_links[j];
        if (lower_bound(path.begin(), path.end(), link1) != path.end()) {
          S[j * N + i] = -1.0;
          sdata.nzn++;
        }
      }
    }

    for (unordered_map<int, unordered_set<int>>::const_iterator it =
             status_primary_path_locs.begin();
         it != status_primary_path_locs.end(); it++) {
      int link = it->first;
      int i = getNindex(link) - K;
      for (unordered_set<int>::const_iterator pit = it->second.begin();
           pit != it->second.end(); pit++) {
        int oindex = paths[*pit].owner;
        for (unordered_set<int>::const_iterator cit =
                 demand_second_path_locs[oindex].begin();
             cit != demand_second_path_locs[oindex].end(); cit++) {
          int slink = paths[*cit].link;
          int j = getNindex(slink) - K;
          S[i * N + j] += 1.0;
          sdata.nzn++;
        }
      }
    }
  }

  void transposeS() {
    double temp = 0.0;
    for (int i = 0; i < N - 1; i++) {
      for (int j = i + 1; j < N; j++) {
        temp = S[i * N + j];
        S[i * N + j] = S[j * N + i];
        S[j * N + i] = temp;
      }
    }
  }

  void computeRHS() {
    fill(X, X + K + N + J, 0.0);
    x_K = X;
    x_N = X + K;
    x_J = X + K + N;

    int nrhs = 1;
    int lda = N;

    int ldb = N;
    int info = 0;
    /*
     * [  I_{K*K}   A            0       ]  [ x_K ] = [ d_K ]
     *bandwidth  envery
     *demand
     * [  B         E            0       ]  [ x_N ] = [ c_N ]
     *capacity of
     *status links
     * [  C         D            I_{J*J} ]  [ x_J ] = [ c_J ]
     *capacity of other
     *links
     * x_N=( B d_K -c_N )/(BA-E)
     *unordered_set S=BA-E, b= B d_K -c_N
     *x_N = b/S
     */
    if (N > 0) {
      fill(b, b + N, 0.0);

      for (int i = 0; i < N; i++) {
        b[i] = -rhs[K + status_links[i]];
      }

      for (int i = 0; i < N; i++) {
        int link = status_links[i];

        const unordered_set<int> &pps = status_primary_path_locs[link];
        for (unordered_set<int>::const_iterator it = pps.begin();
             it != pps.end(); it++) {
          b[i] += rhs[paths[*it].owner];
        }
      }
      copy(S, S + N * N, workS);

      if (KLU == lu_sover) {
        solveLP(workS, N, b);
      } else if (LAPACK == lu_sover) {
        dgesv_(&N, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
               &info);
        sdata.snzn = nonzero(N, workS);

        if (info > 0) {
          assert(false);
          printf(
              "The diagonal element of the "
              "triangular factor of "
              "A,\n");
          printf(
              "U(%i,%i) is zero, so that A is "
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

      memcpy(x_N, b, N * sizeof(double));
    }

    /**
     * x_K=d_K-A x_N
     *
     */
    for (int i = 0; i < K; i++) {
      x_K[i] = rhs[i];
    }

    for (int i = 0; i < K; i++) {
      const unordered_set<int> &pathindices = demand_second_path_locs[i];
      for (unordered_set<int>::const_iterator it = pathindices.begin();
           it != pathindices.end(); it++) {
        x_K[i] -= x_N[getNindex(paths[*it].link) - K];
      }
    }
  }

  /**
   *
   * @brief choose candiate exit base commodity
   *
   * @return
   */
  EXIT_VARIABLE computeLeavingVariable() {
    /**
     *  choose enter base as i which a[ i]>0 and x[ i ]/a[ i ]= min
     * x/a
     *
     */

    int reK = -1;
    double minK = numeric_limits<double>::max();
    double temp;
    int i = 0;
    while (i < K) {
      if (a_K[i] > EPS) {
        reK = i;
        minK = x_K[i] / a_K[i];
        break;
      }
      i++;
    }
    while (i < K) {
      if (a_K[i] > EPS) {
        temp = x_K[i] / a_K[i];
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
    while (i < N) {
      if (a_N[i] > EPS) {
        reN = i;
        minN = x_N[i] / a_N[i];
        break;
      }
      i++;
    }

    while (i < N) {
      if (a_N[i] > EPS) {
        temp = x_N[i] / a_N[i];
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
    while (i < J) {
      if (a_J[i] > EPS) {
        reJ = i;
        minJ = x_J[i] / a_J[i];
        break;
      }
      i++;
    }

    while (i < J) {
      if (a_J[i] > EPS) {
        temp = x_J[i] / a_J[i];
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
      if (minN <= EPS){
        sdata.empty_iterator_num++;
      }
      return re;
    }
    re.id = un_status_links[reJ];
    re.type = OTHER_LINK;
    if (minJ <= EPS){
      sdata.empty_iterator_num++;
    }
    return re;
  }

 public:
  CG(const G &g, const vector<W> &ws, const vector<C> &caps,
     const vector<Demand<C> > &ds)
      : demands(ds),
        orignal_caps(caps),
        thread_num(1),
        lu_sover(KLU) {
    
    orignal_weights.resize( ws.size());    
    vector<int> srcs, snks;
    linkMap.resize(ws.size());
    /**
     *  Reorder link id to improve
     * graph get adjacent nodes speed
     *
     */

    for (size_t v = 0; v < g.getVertex_num(); v++) {
      int degree = g.getOutDegree(v);
      for (int i = 0; i < degree; i++) {
        int link = g.getAdj(v, i);

        int snk;
        linkMap[srcs.size()] = link;
        orignal_weights[srcs.size()]= ws[link];
        
        g.findRhs(link, v, snk);
        srcs.push_back(v);
        snks.push_back(snk);

      }
    }
    graph.initial(srcs, snks);

    
    double max_w = 0;
    for (int i = 0; i < origLink_num; i++) {
      if (max_w < orignal_weights[i]) {
        max_w = orignal_weights[i];
      }
    }
    
    inf_weight = max_w * (graph.getVertex_num() - 1) + 1;

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < 1000; i++) {
      thread_num = omp_get_num_threads();
    }

#endif  // _OPENMP


    info = 0;
    origLink_num = graph.getLink_num();

    status_link_path_loc.resize(origLink_num);

    CZERO = ((C)1e-6);
    EPS = ((double)1e-4);

    K = demands.size();
    demand_second_path_locs.resize(K);

    A = NULL;

    X = NULL;

    ipiv = NULL;

    b = NULL;

    S = NULL;
    workS = NULL;
    S_maxdim = 0;

    update_caps = orignal_caps;
    
    C minNonzero=numeric_limits<C>::max();
    totalB = 0;
    for (typename vector<Demand<C> >::const_iterator it = demands.begin();
         it != demands.end(); it++) {
      if(it->bandwidth>CZERO && it->bandwidth< minNonzero){
        minNonzero=it->bandwidth;
      }
      totalB += it->bandwidth;
    }
    
    if(minNonzero>1){
      for(size_t i=0; i<  update_caps.size(); i++){
        update_caps[i]+=0.0001*(rand()%1000)/1000.0;
      }
      
    }else{
      for(size_t i=0; i<  update_caps.size(); i++){
        update_caps[i]+=minNonzero*(rand()%1000)/1000.0;
      }
    }

    orignal_caps.push_back(inf_weight);
    update_caps.push_back(totalB);
    
  }
  ~CG() {
    if (NULL != A) {
      delete[] A;
      A = NULL;
    }

    if (NULL != X) {
      delete[] X;
      X = NULL;
    }

    if (NULL != ipiv) {
      delete[] ipiv;
      ipiv = NULL;
    }

    if (NULL != b) {
      delete[] b;
      b = NULL;
    }

    if (NULL != S) {
      delete[] S;
      S = NULL;
      delete[] workS;
      workS = NULL;
    }
  }
  void setInfo(const int level) { info = level; }

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
    return (it - un_status_links.begin()) + K + N;
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
    for (int j = 0; j < N + J; j++) {
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

    if (info > 0) {
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

  void initial_solution() {
    paths.resize(K);

    primary_path_loc.resize(K, 0);
    vector<C> temp_p(K);

    orig_success.resize(K, true);
    vector<W> ws(orignal_weights.size(), 1);

#pragma omp parallel for
    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      vector<int> path;
      if (!bidijkstra_shortest_path(graph, ws, src, snk, path,
                                    graph.getLink_num() + 1)) {
        orig_success[i] = false;
      }
    }

    vector<C> temp_cap(orignal_caps);
    vector<bool> success(K, false);

    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      vector<W> ws = orignal_weights;

      for (int j = 0; j < origLink_num; j++) {
        if (temp_cap[j] < bw) {
          ws[j] = inf_weight;
        }
      }

      vector<int> path;
      if (bidijkstra_shortest_path(graph, ws, src, snk, path, inf_weight)) {
        success[i] = true;

        for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
          temp_cap[*it] -= bw;
        }
        sort(path.begin(), path.end());
        paths[i].path = path;
        paths[i].owner = i;
        primary_path_loc[i] = i;
        temp_p[i] = bw;
      }
    }

    
    for (int i = 0; i < K; i++) {
      if (!success[i]) {
        vector<int> path;
        path.push_back(origLink_num);
        paths[i].path = path;
        paths[i].owner = i;
        primary_path_loc[i] = i;
        temp_p[i] = demands[i].bandwidth;
        success[i] = true;
      }
    }
    Inf = inf_weight;

    update_weights = orignal_weights;

    min_commodity_cost.resize(K, 2 * inf_weight);
#pragma omp parallel for
    for (int i = 0; i < K; i++) {
      vector<int> path;
      int src = demands[i].src;
      int snk = demands[i].snk;
      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                                   2 * inf_weight)) {
        min_commodity_cost[i] = path_cost(update_weights, path, (W)0.0);
      }
    }


    
    N = 0;
    J = origLink_num + 1;
    inf_weight = numeric_limits<W>::max() / 3;

    for (int i = 0; i < J; i++) {
      un_status_links.push_back(i);
    }

    rhs.resize(K + origLink_num + 1, (C)0.0);
    for (int i = 0; i < K; i++) {
      rhs[i] = demands[i].bandwidth;
    }

    for (size_t i = 0; i < update_caps.size(); i++) {
      rhs[i + K] = update_caps[i];
    }

    dual_solution.resize(J, 0);

    b = new double[J];
    A = new double[K + J];
    X = new double[K + J];
    fill(X, X + J, 0.0);
    ipiv = new int[J];

    for (int i = 0; i < K; i++) {
      X[i] = temp_p[i];
    }
    x_K = X;
    a_K = A;
  }

  void iteration() {
    while (true) {
      sdata.iterator_num++;
      if (info > 0) {
        if (sdata.iterator_num % 100 == 0) {
          sdata.using_system_time = systemTime() - sdata.start_time;
          C sobj = success_obj();
          std::cout << "using time(s) :" << sdata.using_system_time
                    << std::endl;
          std::cout << sdata.iterator_num << " status link num: " << N
                    << " objvalue: " << computeOBJ()
                    << " success fractional bw: " << sobj
                    << " success rat: " << sobj / (totalB + 0.01)
                    << " the obj gap from opt is: " << sdata.estimee_opt_diff
                    << std::endl;
          std::cout << "nonzero matrix values: " << sdata.nzn
                    << "   SLU non-zeros " << sdata.snzn << std::endl;
          std::cout << sdata.etype << " enter: " << sdata.enter << "    "
                    << sdata.exitt << " exit: " << sdata.exit << std::endl;
        }
      }
      if (sdata.iterator_num > 10000){
        return;
      }

      /**
       *  enter variable choose
       *
       */
      update_edge_left_bandwith();

      ENTER_VARIABLE entering_commodity = chooseEnteringBase();
      if (entering_commodity.id < 0) {
        return;
      }

      EXIT_VARIABLE leaving_base;

      if (PATH_T == entering_commodity.type) {
        /**
         *  exit base  choose
         *
         */
        leaving_base = getLeavingBasebyPath(entering_commodity);
      } else {
        leaving_base = getLeavingBasebyStatusLink(entering_commodity);
      }

      pivot(entering_commodity, leaving_base);

      N = status_links.size();
      J = origLink_num + 1 - N;
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
  ENTER_VARIABLE chooseEnteringBase() {
    ENTER_VARIABLE enter_variable;
    enter_variable.type = PATH_T;
    enter_variable.id = -1;
    /**
     *  check status link dual value
     *
     */
    double max_diff = -EPS;
    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      if (dual_solution[link] < max_diff) {
        max_diff = dual_solution[link];
        enter_variable.id = link;
        enter_variable.type = LINK_T;
        enter_variable.path.clear();
      }
    }

    if (enter_variable.id >= 0) {
      return enter_variable;
    }

    for (size_t i = 0; i < candidate_enter.size(); i++)
      for (vector<int>::iterator it = candidate_enter.begin();
           it != candidate_enter.end(); it++) {
        int id = *it;
        int src = demands[id].src;
        int snk = demands[id].snk;
        vector<int> path;

        W old_cost =
            path_cost(update_weights, paths[primary_path_loc[id]].path, (W)0.0);

        bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                                 inf_weight);
        W new_cost = path_cost(update_weights, path, (W)0.0);

        if (old_cost - new_cost > 10000 * EPS) {
          for (size_t j = 0; j <= i; j++) {
            candidate_enter.erase(candidate_enter.begin());
          }
          enter_variable.id = id;
          enter_variable.type = PATH_T;
          enter_variable.path = path;
          sort(enter_variable.path.begin(), enter_variable.path.end());
          return enter_variable;
        }
      }
    candidate_enter.clear();

    sdata.estimee_opt_diff = 0;
    vector<double> opt_gap(thread_num, 0);
    vector<double> max_diffs(thread_num, EPS);
    vector<double> max_gap(thread_num, EPS);
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
      if (!orig_success[i]) {
        continue;
      }
      int src = demands[i].src;
      int snk = demands[i].snk;
      max_diff = EPS;
      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]].path, (W)0.0);
      if ((old_cost - min_commodity_cost[i]) * demands[i].bandwidth <
          max_diffs[tid]) {
        continue;
      }

      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path[tid],
                                   inf_weight)) {
        W new_cost = path_cost(update_weights, path[tid], (W)0.0);

        if (new_cost > Inf) {
          path[tid].clear();
          path[tid].push_back(origLink_num);
          new_cost = Inf;
        }

        W temp_diff = (old_cost - new_cost);

        if (temp_diff > EPS) {
          if (temp_diff > max_gap[tid]) {
            max_gap[tid] = temp_diff;
          }
          opt_gap[tid] += temp_diff * demands[i].bandwidth;

          temp_diff *= leftBandwith(path[tid]) + EPS;

          if (temp_diff > max_diffs[tid]) {
            if (temp_diff > 10000 * EPS) {
              candidate[tid].push_back(i);
            }

            max_diffs[tid] = temp_diff;

            enter_variables[tid].id = i;
            enter_variables[tid].path = path[tid];
          }
        }
      }
    }

    int chid = 0;
    max_diff = max_diffs[0];
    candidate_enter = candidate[0];
    for (int i = 1; i < thread_num; i++) {
      candidate_enter.insert(candidate_enter.end(), candidate[i].begin(),
                             candidate[i].end());
      if (max_diffs[i] > max_diff) {
        max_diff = max_diffs[i];
        chid = i;
      }
    }

    for (int i = 0; i < thread_num; i++) {
      sdata.estimee_opt_diff += opt_gap[i];
    }

    sort(enter_variables[chid].path.begin(), enter_variables[chid].path.end());

    return enter_variables[chid];
  }

  /**
   * [  I_{K*K}   A            0       ]  [ a_K ] = [ b_K ]
   * [  B         E            0       ]  [ a_N ] = [ b_N ]
   * [  C         D            I_{J*J} ]  [ a_J ] = [ b_J ]
   * a_N=(B b_K-b_N) /( BA-I )
   *
   * [  I_{K*K}   A            0       ]  [ x_K ] = [ d_K ]  bandwidth
   *envery
   *demand
   * [  B         I_{N*N}      0       ]  [ x_N ] = [ c_N ]  capacity of
   *status
   *links
   * [  C         D            I_{J*J} ]  [ x_J ] = [ c_J ]  capacity of
   *other
   *links
   * y_N=( B d_K -c_N )/(BA-I)
   * @param k
   *
   * @return
   */

  EXIT_VARIABLE getLeavingBasebyPath(const ENTER_VARIABLE &enterCommodity) {
    /**
     * A a=d
     *
     */

    fill(A, A + K + N + J, 0.0);
    a_K = A;
    a_N = A + K;
    a_J = A + K + N;

    const vector<int> &path = enterCommodity.path;

    int commodity_primary_path_loc = primary_path_loc[enterCommodity.id];

    const vector<int> &commodity_path = paths[commodity_primary_path_loc].path;

    /**
     *  status links are  empty
     *
     */
    if (0 == N) {
      /**
       * [ I_{K*K}   0       ] [ a_K ]  =[ b_K ]
       * [ B         I_{J*J }] [ a_J ] = [ b_J ]
       *
       *  a_K = b_K
       *  B a_K + a_J = b_J
       *  a_J = b_J -B a_K
       *  a_J= b_J - (b_enterCommodity.id )_J
       */

      a_K[enterCommodity.id] = 1.0;

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        a_J[getJIndex(*it) - N - K] = 1.0;
      }

      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        a_J[getJIndex(*it) - N - K] -= 1.0;
      }

      fill(X, X + K + J, 0);
      x_K = X;
      x_J = X + K;
      /**
       * [ I_{K*K}   0       ] [ x_K ]  =[ rhs_K ]
       * [ B         I_{J*J }] [ x_J ] = [ rhs_J ]
       *
       *  x_K = rhs_K
       *  B rhs_K + x_J = b_J
       *  x_J = b_J -B x_K
       */

      for (int i = 0; i < K; i++) {
        x_K[i] = rhs[i];
      }

      /**
       * y_J = b_J -B y_K
       *
       */

      for (int i = 0; i < J; i++) {
        x_J[i] = edgeLeftBandwith[un_status_links[i]];
      }

    } else {
      /**
       * [  I_{K*K}   A            0       ]  [ a_K ] = [ b_K
       * ]
       * [  B         I_{N*N}      0       ]  [ a_N ] = [ b_N
       * ]
       * [  C         D            I_{J*J} ]  [ a_J ] = [ b_J
       * ]
       * a_N=(B b_K-b_N) /( BA-I )
       */
      int nrhs = 1;
      int lda = N;

      int ldb = N;
      int info = 0;

      /**
       * b=B b_K -b_N
       *
       */

      fill(b, b + N, 0.0);
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
       * a_N=( Bb_K-b_N)/( BA-I )=b/S
       *
       */

      if (KLU == lu_sover) {
        copy(S, S + N * N, workS);
        solveLP(workS, N, b);
      } else if (LAPACK == lu_sover) {
        char c = 'N';
        dgetrs_(&c, &N, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
                &info);

        if (info > 0) {
          assert(false);
          printf(
              "The diagonal element of the "
              "triangular factor of "
              "A,\n");
          printf(
              "U(%i,%i) is zero, so that A is "
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

      memcpy(a_N, b, N * sizeof(double));

      /**
       * a_K=b_K-A a_N
       *
       */

      a_K[enterCommodity.id] = 1.0;
      for (int i = 0; i < K; i++) {
        for (unordered_set<int>::iterator it =
                 demand_second_path_locs[i].begin();
             it != demand_second_path_locs[i].end(); it++) {
          int pindex = *it;
          int link = paths[pindex].link;
          a_K[i] -= a_N[getNindex(link) - K];
        }
      }
      /**
       * a_J=b_J-C a_K-D a_N
       *
       */

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        if (!is_status_link(*it)) a_J[getJIndex(*it) - N - K] = 1.0;
      }

      vector<int> left_AK;
      for (int i = 0; i < K; i++) {
        if (fabs(a_K[i]) > EPS) {
          left_AK.push_back(i);
        }
      }

      for (size_t i = 0; i < un_status_links.size(); i++) {
        for (size_t k = 0; k < left_AK.size(); k++) {
          int pid = left_AK[k];
          int ppid = primary_path_loc[pid];
          if (lower_bound(paths[ppid].path.begin(), paths[ppid].path.end(),
                          un_status_links[i]) != paths[ppid].path.end())
            a_J[i] -= a_K[pid];
        }
      }

      vector<int> left_AN;

      for (int i = 0; i < N; i++) {
        if (fabs(a_N[i]) > EPS) {
          left_AN.push_back(i);
        }
      }

      for (size_t i = 0; i < un_status_links.size(); i++) {
        for (size_t k = 0; k < left_AN.size(); k++) {
          int pid = left_AN[k];
          int ppid = status_link_path_loc[status_links[pid]];
          if (lower_bound(paths[ppid].path.begin(), paths[ppid].path.end(),
                          un_status_links[i]) != paths[ppid].path.end())
            a_J[i] -= a_N[pid];
        }
      }
    }

    return computeLeavingVariable();
  }

  EXIT_VARIABLE getLeavingBasebyStatusLink(const ENTER_VARIABLE &enterLink) {
    /**
     * A a=d
     *
     */

    fill(A, A + K + N + J, 0.0);
    a_K = A;
    a_N = A + K;
    a_J = A + K + N;

    int nrhs = 1;
    int lda = N;

    int ldb = N;
    int info = 0;

    /**
     * b= -b_N
     *
     */
    fill(b, b + N, 0.0);

    b[getNindex(enterLink.id) - K] = -1.0;

    /**
     * a_N=( Bb_K-b_N)/( BA-E )=b/S
     *
     */
    if (KLU == lu_sover) {
      copy(S, S + N * N, workS);
      solveLP(workS, N, b);
    } else if (LAPACK == lu_sover) {
      char c = 'N';
      dgetrs_(&c, &N, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
              &info);
      if (info > 0) {
        assert(false);
        printf(
            "The diagonal element of the triangular "
            "factor of "
            "A,\n");
        printf(
            "U(%i,%i) is zero, so that A is "
            "singular;\n",
            info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
    }

    memcpy(a_N, b, N * sizeof(double));

    /**
     * a_K=-A a_N
     *
     */

    for (int i = 0; i < K; i++) {
      for (unordered_set<int>::iterator it = demand_second_path_locs[i].begin();
           it != demand_second_path_locs[i].end(); it++) {
        int pindex = *it;
        int link = paths[pindex].link;
        assert(link >= 0);
        a_K[i] -= a_N[getNindex(link) - K];
      }
    }

    /**
     * a_J=-C a_K-D a_N
     *
     */

    vector<int> left_AK;
    for (int i = 0; i < K; i++) {
      if (fabs(a_K[i]) > EPS) {
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
          a_J[i] -= a_K[pid];
      }
    }

    vector<int> left_AN;

    for (int i = 0; i < N; i++) {
      if (fabs(a_N[i]) > EPS) {
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
          a_J[i] -= a_N[pid];
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

          assert(!demand_second_path_locs[exit_commodity_id].empty());
          unordered_set<int>::const_iterator it =
              demand_second_path_locs[exit_commodity_id].begin();

          int pid = *it;
          int link = paths[pid].link;

          assert(link >= 0);

          /**
           * from second paths to primary path
           *
           */

          demand_second_path_locs[exit_commodity_id].erase(it);
          primary_path_loc[exit_commodity_id] = pid;
          paths[pid].link = -1;

          demand_second_path_locs[entering_commodity.id].insert(exit_primary_pid);

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
          demand_second_path_locs[paths[pid].owner].erase(pid);
          demand_second_path_locs[entering_commodity.id].insert(pid);
        }
        paths[pid].owner = entering_commodity.id;

      } else {
        // exit is un status link then from un_status
        // link to status link
        if (empty_paths.empty()) {
          Path npath(entering_commodity.path, entering_commodity.id, leaving_base.id);

          paths.push_back(npath);

          status_links.push_back(leaving_base.id);
          status_link_path_loc[leaving_base.id] = paths.size() - 1;

          demand_second_path_locs[entering_commodity.id].insert(paths.size() - 1);

        } else {
          int pid = empty_paths.back();
          empty_paths.pop_back();

          paths[pid].path = entering_commodity.path;
          paths[pid].owner = entering_commodity.id;
          paths[pid].link = leaving_base.id;

          status_links.push_back(leaving_base.id);

          status_link_path_loc[leaving_base.id] = pid;
          demand_second_path_locs[entering_commodity.id].insert(pid);
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
          demand_second_path_locs[leaving_base.id].erase(spid);
        } else {
          assert(!demand_second_path_locs[leaving_base.id].empty());
          unordered_set<int>::const_iterator it =
              demand_second_path_locs[leaving_base.id].begin();

          int pid = *it;
          int link = paths[pid].link;

          assert(link >= 0);

          demand_second_path_locs[leaving_base.id].erase(it);
          primary_path_loc[leaving_base.id] = pid;

          setStatusLink(link, spid);
        }

        addPrimaryPath(exit_commodity_id);
      } else if (STATUS_LINK == leaving_base.type) {
        if (leaving_base.id == entering_commodity.id) {
          empty_paths.push_back(spid);
          demand_second_path_locs[paths[spid].owner].erase(spid);
          paths[spid].owner = -1;
        } else {
          int pid = status_link_path_loc[leaving_base.id];
          empty_paths.push_back(pid);
          demand_second_path_locs[paths[pid].owner].erase(pid);

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

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      allow[paths[pid].owner] += X[i];
      const vector<int> &path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i];
      }
    }

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      allow[paths[pid].owner] += X[i + K];
      const vector<int> &path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i + K];
      }
    }

    for (int i = 0; i < N + J; i++) {
      assert(edgeLeftBandwith[i] > -EPS);
    }
    for (int i = 0; i < K; i++) {
      assert(allow[i] - demands[i].bandwidth > -EPS &&
             allow[i] - demands[i].bandwidth < EPS);
    }

    x_J = X + N + K;
    for (int i = 0; i < J; i++) {
      int link = un_status_links[i];
      x_J[i] = edgeLeftBandwith[link];
      assert(x_J[i] > -EPS);
    }
  }

  C leftBandwith(const vector<int> path) {
    if (path.empty()){
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
   *[ y_K  y_N  0 ]  [  I_{K*K}   A            0       ]   =  [ C_K C_N 0
   *]
   *                 [  B         E            0       ]
   *                 [  C         D            I_{J*J} ]
   *
   * y_K+y_N B=C_K
   * y_K A + y_N =C_N
   *
   * y_N( BA-E ) = C_K A - C_N
   *
   * y_K= C_K -y_N B
   */
  void update_edge_cost() {
    int nrhs = 1;
    int lda = N;

    int ldb = N;
    int info = 0;

    fill(b, b + N, 0.0);
    for (int i = 0; i < N; i++) {
      int linkid = status_links[i];
      int pid = status_link_path_loc[linkid];
      b[i] = -getOrigCost(paths[pid].path);
    }

    fill(A, A + K, 0.0);
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      A[i] = getOrigCost(paths[pid].path);
    }
    for (int i = 0; i < N; i++) {
      int linkid = status_links[i];
      int pid = status_link_path_loc[linkid];
      int oindex = paths[pid].owner;
      b[i] += A[oindex];
    }

    copy(S, S + N * N, workS);
    if (KLU == lu_sover) {
      solveLP(workS, N, b);
    } else if (LAPACK == lu_sover) {
      dgesv_(&N, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb, &info);
      if (info > 0) {
        assert(false);
        printf(
            "The diagonal element of the triangular "
            "factor of "
            "A,\n");
        printf(
            "U(%i,%i) is zero, so that A is "
            "singular;\n",
            info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
    }

    update_weights = orignal_weights;

    fill(dual_solution.begin(), dual_solution.end(), 0.0);

    for (int i = 0; i < N; i++) {
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
              << sobj / (totalB + 0.0) << std::endl;
  }
};

}
}
