
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

#if  (defined(__APPLE__) && defined(__MACH__))

#else
#include <omp.h>

#endif

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include<sstream>
#include<fstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include "System.h"
#include "config.h"
#include "graphalg.hpp"
#include "klu.h"
#include "util.h"


using namespace std;
using milliseconds = std::chrono::duration<double, std::milli>;

namespace raptor {

namespace mcmcf {

struct ENTER_VARIABLE {
  ENTER_BASE_TYPE type;
  int id;
  vector<int> path;
  ENTER_VARIABLE() : type(PATH_T), id(0) {}
};

struct Leaving_base {
  EXIT_BASE_TYPE type;
  int id;
  Leaving_base() : type(DEMAND_T), id(0) {}
};


struct KLUsolver {
  bool first;
  int *Ap;
  int *Ai;
  double *Ax;
  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common Common;
  int nonzeroNum;
  int dim;
  vector<sparseMatrixElem> elements;
  KLUsolver() : first(true), Ap(NULL), Ai(NULL), Ax(NULL) {
    klu_defaults(&Common);
    nonzeroNum = 0;
    dim = 0;
  }
  ~KLUsolver();
  void update(int n);
  bool solve(double *b);
  bool tsolve(double *b);
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

  vector<int> empty_paths;  // the location which  is delete path

  vector<C> dual_solution;  // only link have dual value

  vector<double> min_commodity_cost;

  vector<int> saturate_links;  // contain the index of saturate links

  vector<int> un_saturate_links;

  vector<int> un_saturate_link_ids;// the id of un saturate link
  

  vector<unordered_set<int>>
      demand_secondary_path_locs;  // belong to fixed demands' second paths

  unordered_map<int, unordered_set<int>>
      saturate_primary_path_locs;  // primary paths which
                                   // corross the saturate link

  vector<int> primary_path_loc;  // every demands have a primary path

  vector<int> saturate_link_path_loc;  //  the path corresponding saturate link

  Statistics_data sdata;

  vector<int> candidate_enter;

  KLUsolver klusolver;

  vector<C> rhs;

  double *X;
  double *Lambda;

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
  vector<pair<int, int> >  saturateLinkAndMatrix;

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

      if (paths[pid].back() < origLink_num) {
        double p_cost = path_cost(orignal_weights, paths[pid], 0.0);

        re += p_cost * x_K[i];
      }
    }

    for (int i = 0; i < S; i++) {
      int link = saturate_links[i];
      int pid = saturate_link_path_loc[link];
      if (paths[pid].back() < origLink_num) {
        double p_cost = path_cost(orignal_weights, paths[pid], 0.0);

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
   * column primary
   *
   */
  void computeS() {
    if (0 == S) {
      return;
    }
    if (KLU == para.solver) {
      vector<sparseMatrixElem> elements;

      /**
       * -D
       */
      for (int i = 0; i < S; i++) {
        int link = saturate_links[i];
        int pid = saturate_link_path_loc[link];
        const vector<int> &path = paths[pid].path;
        for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
             lit++) {
          int j = binfind(saturate_links.begin(), saturate_links.end(), *lit);
          if (j > -1) {
            sparseMatrixElem elem;
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
      for (unordered_map<int, unordered_set<int>>::const_iterator it =
               saturate_primary_path_locs.begin();
           it != saturate_primary_path_locs.end(); it++) {
        int link = it->first;
        int i = getSindex(link) - K;
        for (unordered_set<int>::const_iterator pit = it->second.begin();
             pit != it->second.end(); pit++) {
          int oindex = paths[*pit].owner;
          for (unordered_set<int>::const_iterator cit =
                   demand_secondary_path_locs[oindex].begin();
               cit != demand_secondary_path_locs[oindex].end(); cit++) {
            int slink = paths[*cit].link;
            int j = getSindex(slink) - K;

            sparseMatrixElem elem;
            elem.column = j;
            elem.row = i;
            elem.value = 1.0;
            elements.push_back(elem);
          }
        }
      }

      klusolver.elements = elements;
      klusolver.update(S);
      sdata.nzn = klusolver.nonzeroNum;

    } else {
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
          int j = binfind(saturate_links.begin(), saturate_links.end(), *lit);
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
      for (unordered_map<int, unordered_set<int>>::const_iterator it =
               saturate_primary_path_locs.begin();
           it != saturate_primary_path_locs.end(); it++) {
        int link = it->first;
        int i = getSindex(link) - K;
        for (unordered_set<int>::const_iterator pit = it->second.begin();
             pit != it->second.end(); pit++) {
          int oindex = paths[*pit].owner;
          for (unordered_set<int>::const_iterator cit =
                   demand_secondary_path_locs[oindex].begin();
               cit != demand_secondary_path_locs[oindex].end(); cit++) {
            int slink = paths[*cit].link;
            int j = getSindex(slink) - K;
            SM[i * S + j] += 1.0;
            sdata.nzn++;
          }
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
    fill(X, X + K + S + N, 0.0);
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
      for (int i = 0; i < S; i++) {
        b[i] = -rhs[K + saturate_links[i]];
      }
      /**
       * C rhs_K
       *
       */

      for (int i = 0; i < S; i++) {
        int link = saturate_links[i];

        const unordered_set<int> &pps = saturate_primary_path_locs[link];
        for (unordered_set<int>::const_iterator it = pps.begin();
             it != pps.end(); it++) {
          b[i] += rhs[paths[*it].owner];
        }
      }
      double t=0;
      if (KLU == para.solver) {

        callTime(t,klusolver.solve(b));

      } else if (LAPACK == para.solver) {
        copy(SM, SM + S * S, workS);
        callTime(t,dgesv_(&S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
                          &info));
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
      sdata.lpsolvertime+=t;
      memcpy(x_S, b, S * sizeof(double));
      for (int i = 0; i < S; i++) {
        int link = saturate_links[i];
        int spid = saturate_link_path_loc[link];
        assert(b[i] <= rhs[paths[spid].owner]);
      }
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
        x_K[i] -= x_S[getSindex(paths[*it].link) - K];
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
   * @brief choose candiate leaving base
   *
   * @return
   */
  Leaving_base computeLeavingVariable() {
    /**
     *  choose entering base as i=argmin_{ lambda[ i]>0} x[ i ]/lambda[ i ]
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

    int reS = -1;
    double minS = numeric_limits<double>::max();
    i = 0;
    while (i < S) {
      if (lambda_S[i] > EPS) {
        reS = i;
        minS = x_S[i] / lambda_S[i];
        break;
      }
      i++;
    }

    while (i < S) {
      if (lambda_S[i] > EPS) {
        temp = x_S[i] / lambda_S[i];
        if (temp < minS) {
          reS = i;
          minS = temp;
        }
      }
      i++;
    }

    int reN = -1;
    double minN = numeric_limits<double>::max();

    i = 0;
    while (i < N) {
      if (lambda_N[i] > EPS) {
        reN = i;
        minN = x_N[i] / lambda_N[i];
        break;
      }
      i++;
    }

    while (i < N) {
      if (lambda_N[i] > EPS) {
        temp = x_N[i] / lambda_N[i];
        if (temp < minN) {
          reN = i;
          minN = temp;
        }
      }
      i++;
    }

    Leaving_base re;

    if (minK <= minS && minK <= minN) {
      re.id = reK;
      re.type = DEMAND_T;
      if (minK <= EPS/100.0){
        sdata.empty_iterator_num++;
      }
      return re;
    }
    if (minS <= minN) {
      re.id = saturate_links[reS];
      re.type = STATUS_LINK;
      if (minS <= EPS/100.0) {
        sdata.empty_iterator_num++;
      }
      return re;
    }
    re.id = un_saturate_links[reN];
    re.type = OTHER_LINK;
    if (minN <= EPS/100.0) {
      sdata.empty_iterator_num++;
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

    /**
     *  Reorder links' id to improve
     * graph get adjacent nodes speed
     *
     */

    for (size_t v = 0; v < g.getVertex_num(); v++) {
      int degree = g.getOutDegree(v);
      for (int i = 0; i < degree; i++) {
        int link = g.getAdj(v, i);
        int id=srcs.size();
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
    if(para.isSetpenaltyPrice){
      
      inf_weight=para.penaltyPriceForFailDemand;
      
    }else{
      inf_weight = max_w * sqrt(graph.getVertex_num()) + 1;
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

    if(para.isSetDisturbed){
      double domain=pow(10,para.disturbedDigit)+1;
      double multiRat=1.0/domain;
      
      if (minNonzero > 1) {
        for (size_t i = 0; i < update_caps.size(); i++) {
          update_caps[i] += multiRat * (rand() /( RAND_MAX+0.0));
        }

      } else {
        for (size_t i = 0; i < update_caps.size(); i++) {
          update_caps[i] += multiRat*minNonzero * (rand() /( RAND_MAX+0.0));
        }
      }
    }
    /**
     *add a  dummy link which can setup all demands
     * 
     */

    orignal_weights.push_back(inf_weight);
    update_caps.push_back(TotalB+1);
    inf_weight=getInf(0.0);
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

    if (NULL != SM) {
      delete[] SM;
      SM = NULL;
      delete[] workS;
      workS = NULL;
    }
  }
  void setPara(solverPara &p){
    para=p;
  }
      
  void setInfo(const int level) { para.info = level; }

  void setLUSOLVER(LU_SOLVER s) { para.solver = s; }



  void setStatusLink(const int link, const int pid) {
    paths[pid].link = link;
    saturate_link_path_loc[link] = pid;
  }

  void addPrimaryPath(const int commodityID) {
    int pid = primary_path_loc[commodityID];
    const vector<int> &new_path = paths[pid].path;
    for (vector<int>::const_iterator it = new_path.begin();
         it != new_path.end(); it++) {
      if (binfind(saturate_links.begin(), saturate_links.end(), *it) > -1) {
        saturate_primary_path_locs[*it].insert(pid);
      }
    }
  }

  void deletePrimaryPath(const int commodityID) {
    int primary_pid = primary_path_loc[commodityID];
    const vector<int> &orig_path = paths[primary_pid].path;
    for (vector<int>::const_iterator it = orig_path.begin();
         it != orig_path.end(); it++) {
      if (binfind(saturate_links.begin(), saturate_links.end(), *it) > -1) {
        saturate_primary_path_locs[*it].erase(primary_pid);
      }
    }
  }

  void addStatusLink(const int link) {
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      if (find(paths[pid].path.begin(), paths[pid].path.end(), link) !=
          paths[pid].path.end()) {
        saturate_primary_path_locs[link].insert(pid);
      }
    }
  }
  void deleteStatusLink(const int link) {
    vector<int>::iterator it =
        lower_bound(saturate_links.begin(), saturate_links.end(), link);
#if DEBUG
    assert(*it == link && it != saturate_links.end());
#endif

    saturate_links.erase(it);

    int spid = saturate_link_path_loc[link];
    paths[spid].link = -1;

    saturate_primary_path_locs.erase(link);
  }

  /**
   *
   *
   * @param i the id of saturate link
   *
   * @return the index of the saturate link in simplex matrix
   */
  int getSindex(int i) const {
#ifdef DEBUG
    assert(lower_bound(saturate_links.begin(), saturate_links.end(), i) !=
           saturate_links.end());
#endif

    int j = binfind(saturate_links.begin(), saturate_links.end(), i);
    return j + K;
  }

  double getOrigCost(const vector<int> &path) const {
    return path_cost(orignal_weights, path, 0.0);
  }

  void computIndexofLinks() {

    un_saturate_links.clear();

    fill(un_saturate_link_ids.begin(), un_saturate_link_ids.end(), -1);
    size_t i = 0;
    for (int j = 0; j < S + N; j++) {
      if (i >= saturate_links.size()) {
        un_saturate_link_ids[j]=un_saturate_links.size();
        un_saturate_links.push_back(j);
      } else {
        if (j < saturate_links[i]) {
          un_saturate_link_ids[j]=un_saturate_links.size();
          un_saturate_links.push_back(j);
        } else {
          assert(j==saturate_links[i]);
          i++;
        }
      }
    }
  }

  bool solve() {
    sdata.start_time = systemTime();
    initial_solution();

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

    vector<pair<C, int> > sortDemands;
    for (int i = 0; i < K; i++) {
      pair<C, int> temp =make_pair(demands[i].bandwidth, i );
      sortDemands.push_back(temp);
    }

    sort(sortDemands.rbegin(), sortDemands.rend());

    /**
     * Direcltly setup demand one by one
     *
     */

    for (int k = 0; k < K; k++) {
      
      int i=sortDemands[k].second;
      
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      vector<W> ws( orignal_weights.size(), 1);

      for (int j = 0; j < origLink_num; j++) {
        if (temp_cap[j] < bw) {
          ws[j] = graph.getVertex_num();
        }else{
          ws[j]=bw/temp_cap[j];
        }
        
      }

      vector<int> path;
      if (bidijkstra_shortest_path(graph, ws, src, snk, path, graph.getVertex_num())) {
        success[i] = true;

        for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
          temp_cap[*it] -= bw;
        }

        paths[i].path = path;
        paths[i].owner = i;
        primary_path_loc[i] = i;
      }
    }

    for (int i = 0; i < K; i++) {
      if (!success[i]) {
        vector<int> path;
        path.push_back(origLink_num);
        paths[i].path = path;
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
    
    for (int i = 0; i < N; i++) {
      un_saturate_links.push_back(i);
      un_saturate_link_ids.push_back(i);
    }

    rhs.resize(K + origLink_num + 1, (C)0.0);
    for (int i = 0; i < K; i++) {
      rhs[i] = demands[i].bandwidth;
    }

    for (size_t i = 0; i < update_caps.size(); i++) {
      rhs[i + K] = update_caps[i];
    }

    dual_solution.resize(N + 1, 0);

  }

  bool iteration() {
    while (true) {
      sdata.iterator_num++;
      computeRHS();
      double OBJ=computeOBJ();
      sdata.totalStaturateLink+=S;
      sdata.totalNonzero+=sdata.nzn;
      if(para.info>1){
        saturateLinkAndMatrix.push_back(make_pair(S,sdata.nzn ));
      }

      if(sdata.bestUpperobj<OBJ-sdata.estimee_opt_diff){
        sdata.bestUpperobj=OBJ-sdata.estimee_opt_diff;
      }
      if (sdata.iterator_num % para.perIterationPrint == 0)
        cout<<paths.size()<<endl;
      if (para.info > 2) {
        if (sdata.iterator_num % para.perIterationPrint == 0) {
          sdata.using_system_time = systemTime() - sdata.start_time;
          C sobj = success_obj();
          
          std::cout << fixed;
          std::cout << "============================================ " << endl;
          std::cout << "iteration: " << sdata.iterator_num << endl;

          std::cout << "using time(s) :" << sdata.using_system_time
                    << std::endl;
          std::cout << "empty iteration num: " << sdata.empty_iterator_num
                    << endl;
          std::cout << "saturate link num: " << S << endl;
          std::cout << "objvalue: " << OBJ << endl;
          std::cout << "success fractional bw: " << sobj
                    << ", success rat: " << sobj / (TotalB + 0.01)<< std::endl;
          std::cout << "The obj gap from opt less or equal than: " << sdata.estimee_opt_diff<<"("<<100*sdata.estimee_opt_diff/sdata.bestUpperobj<<"%)"
                    << std::endl;

          std::cout << "The number of nonzeon element in matrix (CB-D) : "
                    << sdata.nzn << std::endl;

          std::cout << "Last entering: " << sdata.enter
                    << ", last leaving: " << sdata.exit << std::endl;
        }
      }
      if (sdata.iterator_num > para.maxIterationNum) {
        return false;
      }

      update_edge_left_bandwith();
      /**
       *  entering variable choose
       *
       */

      double t=0;
      callTime(t, ENTER_VARIABLE entering_commodity = chooseEnteringVariable());
      sdata.shortestpathtime+=t;
      
      if (entering_commodity.id < 0) {
        return true;
      }

      Leaving_base leaving_base;
      /**
       *  leaving base  choose
       *
       */
      if (PATH_T == entering_commodity.type) {
        leaving_base = getLeavingBasebyPath(entering_commodity);
      } else {
        leaving_base = getLeavingBasebyStatusLink(entering_commodity);
      }

      pivot(entering_commodity, leaving_base);

      S = saturate_links.size();
      assert(empty_paths.size()+K+S==paths.size());
      N = origLink_num + 1 - S;
      computIndexofLinks();

      computeS();

      update_edge_cost();

      /**
       * column primary
       *
       */
      if (KLU != para.solver) {
        transposeS();
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
  ENTER_VARIABLE chooseEnteringVariable() {
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
      if (dual_solution[link] < max_diff) {
        max_diff = dual_solution[link];
        enter_variable.id = link;
      }
    }

    /**
     * If there is a dual value of  saturate link is negative then this link is
     * a entering variable
     */

    if (enter_variable.id >= 0) {
      return enter_variable;
    }

    /**
     * fast way to choose entering variable as a path
     */
    for (size_t i = 0; i < candidate_enter.size(); i++) {
      int id = candidate_enter[i];
      int src = demands[id].src;
      int snk = demands[id].snk;
      vector<int> path;

      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[id]].path, (W)0.0);

      if (bidijkstra_shortest_path(graph, update_weights, src, snk, path,
                                   inf_weight)) {
        W new_cost = path_cost(update_weights, path, (W)0.0);
        W diff = old_cost - new_cost;

        if (diff > 10000 * EPS && diff > sdata.objSpeed&& leftBandwith(path)>(demands[id].bandwidth/10.0)) {
          candidate_enter.erase(candidate_enter.begin(),
                                candidate_enter.begin() + i + 1);
          enter_variable.id = id;
          enter_variable.type = PATH_T;
          enter_variable.path = path;

          sdata.objSpeed = (para.objSpeedUpdateRat) * sdata.objSpeed +
                           (1 - (para.objSpeedUpdateRat)) * diff;
          return enter_variable;
        }
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
      if (!orig_success[i]) {
        continue;
      }
      int src = demands[i].src;
      int snk = demands[i].snk;
      max_diff = EPS;
      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]].path, (W)0.0);
      /**
       * If the possible biggest reduce objective value is less than exist found
       * one then continue
       *
       */

      if ((old_cost - min_commodity_cost[i]) * demands[i].bandwidth <
          max_gaps[tid]) {
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

        if (temp_diff > EPS) {
          if (temp_diff > max_diffs[tid]) {
            max_diffs[tid] = temp_diff;
          }
          opt_gap[tid] += temp_diff * demands[i].bandwidth;

          temp_diff;  //*= leftBandwith(path[tid]) + EPS;

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

    sdata.objSpeed = (para.objSpeedUpdateRat) * sdata.objSpeed +
                     (1 - (para.objSpeedUpdateRat)) * max_diffs[chid];
    return enter_variables[chid];
  }

  /**
   * data_N -=  -H data_K-F data_S
   *
   */
  void computeLeftN(double *data) {
    double *data_K = data;
    double *data_S = data + K;
    double *data_N = data + K + S;

    vector<int> left_dataK;
    for (int i = 0; i < K; i++) {
      if (fabs(data_K[i]) > EPS) {
        left_dataK.push_back(i);
      }
    }

    for (size_t k = 0; k < left_dataK.size(); k++) {
      int pid = left_dataK[k];
      int ppid = primary_path_loc[pid];
      for (vector<int>::const_iterator lit = paths[ppid].path.begin();
           lit != paths[ppid].path.end(); lit++) {
        int i =un_saturate_link_ids[*lit];

        if (i > -1) {
          data_N[i] -= data_K[pid];
        }
      }
    }

    vector<int> left_dataS;

    for (int i = 0; i < S; i++) {
      if (fabs(data_S[i]) > EPS) {
        left_dataS.push_back(i);
      }
    }

    for (size_t k = 0; k < left_dataS.size(); k++) {
      int lid = left_dataS[k];
      int ppid = saturate_link_path_loc[saturate_links[lid]];
      for (vector<int>::const_iterator lit = paths[ppid].path.begin();
           lit != paths[ppid].path.end(); lit++) {
        int i =un_saturate_link_ids[*lit];

        if (i > -1) {
          data_N[i] -= data_S[lid];
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
   * @return leaving base
   */

  Leaving_base getLeavingBasebyPath(const ENTER_VARIABLE &enterCommodity) {
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

      fill(b, b + S, 0.0);
      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        int i = binfind(saturate_links.begin(), saturate_links.end(), *it);

        if (i > -1) {
          b[i] = 1.0;
        }
      }

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        int i = binfind(saturate_links.begin(), saturate_links.end(), *it);

        if (i > -1) {
          b[i] -= 1.0;
        }
      }
      /**
       * lambda_S=( C beta_K-beta_S)/( CB - D )=b/SM
       *
       */
      double t=0;
      if (KLU == para.solver) {
        callTime(t, klusolver.solve(b));
      } else if (LAPACK == para.solver) {
        copy(SM, SM + S * S, workS);
        callTime(t,dgesv_(&S, &nrhs, workS, &lda, ipiv, b, &ldb, &info));

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
      sdata.lpsolvertime+=t;
      memcpy(lambda_S, b, S * sizeof(double));

      /**
       * lambda_K=beta_K-B lambda_S
       *
       */

      lambda_K[enterCommodity.id] = 1.0;
      for (int i = 0; i < K; i++) {
        for (unordered_set<int>::iterator it =
                 demand_secondary_path_locs[i].begin();
             it != demand_secondary_path_locs[i].end(); it++) {
          int pindex = *it;
          int link = paths[pindex].link;
          lambda_K[i] -= lambda_S[getSindex(link) - K];
        }
      }

      /**
       * lambda_N=beta_N- H lambda_K-F lambda_S
       *
       */

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        int i =un_saturate_link_ids[*it];

        if (i > -1) {
          lambda_N[i] = 1.0;
        }
      }

      computeLeftN(Lambda);
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
  Leaving_base getLeavingBasebyStatusLink(const ENTER_VARIABLE &enterLink) {
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

    b[getSindex(enterLink.id) - K] = -1.0;

    /**
     * lambda_S=( C beta_K-beta_S)/( CB - D )=b/SM
     *
     */
    double t=0;
    if (KLU == para.solver) {
      callTime(t,klusolver.solve(b));

    } else if (LAPACK == para.solver) {
      copy(SM, SM + S * S, workS);
      callTime(t,dgesv_(&S, &nrhs, workS, &lda, ipiv, b, &ldb, &info));
      // char c = 'S';
      // dgetrs_(&c, &S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb,
      //         &info);
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
    sdata.lpsolvertime+=t;
    memcpy(lambda_S, b, S * sizeof(double));

    /**
     * lambda_K=-B lambda_S
     *
     */

    for (int i = 0; i < K; i++) {
      for (unordered_set<int>::iterator it =
               demand_secondary_path_locs[i].begin();
           it != demand_secondary_path_locs[i].end(); it++) {
        int pindex = *it;
        int link = paths[pindex].link;
        assert(link >= 0);
        lambda_K[i] -= lambda_S[getSindex(link) - K];
      }
    }

    /**
     * lambda_N =  -H lambda_K-F lambda_S
     *
     */
    computeLeftN(Lambda);

    return computeLeavingVariable();
  }

  void pivot(ENTER_VARIABLE &entering_commodity, Leaving_base &leaving_base) {
    sdata.enter = entering_commodity.id;
    sdata.exit = leaving_base.id;
    if (PATH_T == entering_commodity.type) {
      if (DEMAND_T == leaving_base.type) {
        int exit_commodity_id = leaving_base.id;
        int exit_primary_pid = primary_path_loc[exit_commodity_id];

        /**
         * leaving primary will been deleted from base
         * matrix
         *
         */
        deletePrimaryPath(exit_commodity_id);

        /**
         * when entering commodity and leaving commodity are
         *same then replace the
         *commodity primary path with  entering path
         *
         */

        paths[exit_primary_pid].path = entering_commodity.path;
        paths[exit_primary_pid].owner = entering_commodity.id;

        if (entering_commodity.id != leaving_base.id) {
          /**
           * when entering commodity and leaving
           *commodity are diff then replace the
           *commodity primary path with the second
           *path of the commodity which
           * crossponding saturate link cross entering
           *path and make entering path as
           *path coreesponding saturate link
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
        int pid = saturate_link_path_loc[leaving_base.id];
        paths[pid].path = entering_commodity.path;

        /**
         * when the owner of the leaving path is not owner
         * of entering path
         *
         */

        if (entering_commodity.id != paths[pid].owner) {
          demand_secondary_path_locs[paths[pid].owner].erase(pid);
          demand_secondary_path_locs[entering_commodity.id].insert(pid);
        }
        paths[pid].owner = entering_commodity.id;

      } else {
        // leaving is un saturate link then from un_saturate
        // link to saturate link
        if (empty_paths.empty()) {
          Path npath(entering_commodity.path, entering_commodity.id,
                     leaving_base.id);

          paths.push_back(npath);

          saturate_links.push_back(leaving_base.id);
          stable_sort(saturate_links.begin(), saturate_links.end());
          saturate_link_path_loc[leaving_base.id] = paths.size() - 1;

          demand_secondary_path_locs[entering_commodity.id].insert(
              paths.size() - 1);

        } else {
          int pid = empty_paths.back();
          empty_paths.pop_back();

          paths[pid].path = entering_commodity.path;
          paths[pid].owner = entering_commodity.id;
          paths[pid].link = leaving_base.id;

          saturate_links.push_back(leaving_base.id);
          stable_sort(saturate_links.begin(), saturate_links.end());
          saturate_link_path_loc[leaving_base.id] = pid;
          demand_secondary_path_locs[entering_commodity.id].insert(pid);
        }

        int link = leaving_base.id;
        addStatusLink(link);
      }

    } else {
      /**
       * entering a saturate link
       *
       */
      int enter_status_link = entering_commodity.id;

      int spid = saturate_link_path_loc[enter_status_link];
      deleteStatusLink(enter_status_link);

      if (DEMAND_T == leaving_base.type) {
        int exit_commodity_id = leaving_base.id;

        deletePrimaryPath(exit_commodity_id);
        empty_paths.push_back(primary_path_loc[leaving_base.id]);
        paths[primary_path_loc[leaving_base.id]].path.clear();
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
          paths[spid].path.clear();
        } else {
          int pid = saturate_link_path_loc[leaving_base.id];
          empty_paths.push_back(pid);
          paths[pid].path.clear();
          demand_secondary_path_locs[paths[pid].owner].erase(pid);

          paths[pid].owner = -1;
          paths[pid].link = -1;

          setStatusLink(leaving_base.id, spid);
        }

      } else {
        saturate_links.push_back(leaving_base.id);
        stable_sort(saturate_links.begin(), saturate_links.end());
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

    fill(b, b + S, 0.0);
    for (int i = 0; i < S; i++) {
      int linkid = saturate_links[i];
      int pid = saturate_link_path_loc[linkid];
      b[i] = getOrigCost(paths[pid].path);
    }

    for (int i = 0; i < S; i++) {
      int linkid = saturate_links[i];
      int pid = saturate_link_path_loc[linkid];
      int oindex = paths[pid].owner;
      int ppid = primary_path_loc[oindex];
      b[i] -= getOrigCost(paths[ppid].path);
    }
    double t=0;
    if (KLU == para.solver) {
      callTime(t,klusolver.tsolve(b));
    } else if (LAPACK == para.solver) {
      copy(SM, SM + S * S, workS);
      callTime(t,dgesv_(&S, &nrhs, (double *)workS, &lda, ipiv, (double *)b, &ldb, &info));
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
    sdata.lpsolvertime+=t;
    update_weights = orignal_weights;

    fill(dual_solution.begin(), dual_solution.end(), 0.0);

    for (int i = 0; i < S; i++) {
      dual_solution[saturate_links[i]] = b[i];
      update_weights[saturate_links[i]] += b[i];
    }
  }

  void printResult() {
    if(para.info>1){
      sdata.using_system_time = systemTime() - sdata.start_time;
      C sobj = success_obj();
      std::cout<<systemTime() - sdata.start_time<<","<<sobj / (TotalB + 0.0)<<","<<sobj<<","<< sdata.iterator_num<<","
               <<sdata.empty_iterator_num <<","<<sdata.shortestpathtime/1000.0<<","<<sdata.lpsolvertime<<std::endl;
            
    }else{
      
      sdata.using_system_time = systemTime() - sdata.start_time;
      std::cout << "======================================"<<std::endl;
      std::cout << "using time(s) :" << sdata.using_system_time << std::endl;
      std::cout << "computing shortest path use time(ms) :" << sdata.shortestpathtime << std::endl;
      std::cout << "solving linear equation solve use time(ms) :" << sdata.lpsolvertime << std::endl;
      std::cout << "iteration time: " << sdata.iterator_num << std::endl;
      std::cout << "empty iteration tiem: " << sdata.empty_iterator_num
                << std::endl;
      std::cout<<"Total nonzero element: Total saturate link "<<sdata.totalNonzero/(sdata.totalStaturateLink+0.01)<<std::endl;

      C sobj = success_obj();

      std::cout << "success fractional bandwidth: " << sobj << std::endl;
      std::cout << "success fractional bandwidth rat in total demand: "
                << sobj / (TotalB + 0.0) << std::endl;
    }
    if(para.info>1){
      std::fstream fs;
      stringstream sstr;
      if(KLU==para.solver){
        sstr<<"klusaturate"<<graph.getVertex_num()<<"_"<<graph.getLink_num()<< ".csv";
      }else{
        sstr<<"blassaturate"<<graph.getVertex_num()<<"_"<<graph.getLink_num()<< ".csv";
      }

      fs.open (sstr.str(),  std::fstream::out|std::fstream::trunc);
      for(vector<pair<int, int> >::iterator it=saturateLinkAndMatrix.begin(); it!= saturateLinkAndMatrix.end(); it++){
        fs<<it->first<<","<<it->second<<endl;
      }
      fs.close();
      
    }
  }
};
}
}
