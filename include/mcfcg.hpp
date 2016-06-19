
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
#include <limits>
#include <set>
#include <map>
#include <cstring>
#include <cassert>
#include <limits>
#include <algorithm>
#include <iostream>
#include<omp.h>
// #include "klu.h"
#include "graphalg.hpp"
using namespace fast_graph;
using namespace std;

/* DGESV prototype */
extern "C" {
  void dgesv_(int* N, int* nrhs, double* a, int* lda, int* ipiv, double* b,
              int* ldb, int* info);
  
  void sgesv_(int* N, int* nrhs, float* a, int* lda, int* ipiv, float* b,
                int* ldb, int* info);

  // LU decomoposition of a general matrix
  void dgetrf_( int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  //// generate inverse of a matrix given its LU decomposition
  //void dgetri_( int* N, double* A, int* lda, int* IPIV, double* WORK, int*    lwork, int* INFO);
  void dgetrs_( char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);


    // LU decomoposition of a general matrix
  void sgetrf_( int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);
  //// generate inverse of a matrix given its LU decomposition
  //void dgetri_( int* N, double* A, int* lda, int* IPIV, double* WORK, int*    lwork, int* INFO);
  void sgetrs_( char* C, int* N, int* NRHS, float* A, int* LDA, int* IPIV, float* B, int* LDB, int* INFO);

  
}

namespace mcmcf {
enum ENTER_BASE_TYPE {
  PATH_T = 0,
  LINK_T = 1

};

struct ENTER_VARIABLE {
  ENTER_BASE_TYPE type;
  int id;
  vector<int> path;
  ENTER_VARIABLE() : type(PATH_T), id(0) {}
};

enum EXIT_BASE_TYPE {
  DEMAND_T = 0,
  STATUS_LINK = 1,
  OTHER_LINK = 2

};
struct EXIT_VARIABLE {
  EXIT_BASE_TYPE type;
  int id;
  EXIT_VARIABLE() : type(DEMAND_T), id(0) {}
};

template <typename G, typename W, typename C>
class CG {
 public:
  struct Demand {
    int src, snk;
    C bandwidth;
    Demand() : src(0), snk(0), bandwidth(0) {}
  };

  struct Statistics_data {
    int iterator_num;
    int empty_iterator_num;
    W estimee_opt_diff;
    Statistics_data() : iterator_num(0), empty_iterator_num(0),estimee_opt_diff( numeric_limits<W>::max(  ) ) {}
  };

  struct Path {
    vector<int> path;
    int owner;
    int link;

    Path() : owner(-1), link(-1) {}
    Path(const vector<int>& p, const int o, const int l)
        : path(p), owner(o), link(l) {}
  };

 private:
  G graph;
  int origLink_num;
  G newGraph;
  W inf_weight;

  vector<Demand> demands;

  vector<W> orignal_weights;
  vector<C> orignal_caps;

  vector<W> update_weights;
  vector<C> update_caps;
  vector<C> edgeLeftBandwith;

  C totalB;

  vector<Path> paths;  // all the paths save in this vector
  vector<int> empty_paths;  // the location which  is delete path
  // vector<int> owner; // every path has over demand
  // vector<int> link_of_path;

  vector<C> dual_solution;  // only link have dual value

  vector<int> status_links;  // contain the index of status links

  vector<int> un_status_links;

  vector<set<int> >
  demand_second_path_locs;  // belong to fixed demands' second paths

  map<int, set<int> >
  status_primary_path_locs;  // primary paths which corross the status link

  vector<int> primary_path_loc;  // every demands have a primary path

  map<int, int> status_link_path_loc;  //  the path corresponding status link

  Statistics_data sdata;

  vector<C> rhs;
  W* A;
  W* X;

  int* ipiv;
  W* b;

  W* S;
  W* workS;
  int S_maxdim;

  W* a_K;
  W* a_N;
  W* a_J;

  W* x_K;
  W* x_N;
  W* x_J;

  W* y_K;
  W* y_N;
  W* y_J;

  C CZERO;
  int K, N, J;
  W EPS;

  int info;

  int thread_num;

  void allocateS(const int N) {
    if (S_maxdim >= N) return;
    if (NULL != S) {
      delete[] S;
      S = NULL;
      delete[] workS;
    }
    S_maxdim = 1.3 * N + 1;
    S = new W[S_maxdim * S_maxdim];
    workS = new W[S_maxdim * S_maxdim];
  }
  C success_obj() {
    C re = 0;
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      if (paths[pid].path.back() < origLink_num) re += x_K[i];
    }

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];

      if (paths[pid].path.back() < origLink_num) re += x_N[i];
    }
    return re;
  }
  W computeOBJ() {
    W re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      W p_cost = 0;
      for (vector<int>::const_iterator it = paths[pid].path.begin();
           it != paths[pid].path.end(); it++) {
        p_cost += orignal_weights[*it];
      }
      re += p_cost * x_K[i];
    }

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];

      W p_cost = 0;
      for (vector<int>::const_iterator it = paths[pid].path.begin();
           it != paths[pid].path.end(); it++) {
        p_cost += orignal_weights[*it];
      }
      re += p_cost * x_N[i];
    }
    return re;
  }
  W computeObj() {
    W re = 0;

    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];

      if (paths[pid].back() < origLink_num) {
        W p_cost = 0;
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
        W p_cost = 0;
        for (vector<int>::const_iterator it = paths[pid].begin();
             it != paths[pid].end(); it++) {
          p_cost += orignal_weights[*it];
        }
        re += p_cost * x_N[i];
      }
    }
    return re;
  }

  /**
   * row primary
   *
   */
  void computeS() {
    if (0 == N) return;
    /**
     * S=BA-E
     *
     */
    if (N > S_maxdim) allocateS(N);
    fill(S, S + N * N, 0.0);

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      const vector<int>& path = paths[pid].path;
      for (int j = 0; j < N; j++) {
        int link1 = status_links[j];
        if (find(path.begin(), path.end(), link1) != path.end()) {
          S[j * N + i] = -1.0;
        }
      }
    }

    for (map<int, set<int> >::const_iterator it =
             status_primary_path_locs.begin();
         it != status_primary_path_locs.end(); it++) {
      int link = it->first;
      int i = getNindex(link) - K;
      for (set<int>::const_iterator pit = it->second.begin();
           pit != it->second.end(); pit++) {
        int oindex = paths[*pit].owner;
        for (set<int>::const_iterator cit =
                 demand_second_path_locs[oindex].begin();
             cit != demand_second_path_locs[oindex].end(); cit++) {
          int slink = paths[*cit].link;
          int j = getNindex(slink) - K;
          S[i * N + j] += 1.0;
        }
      }
    }
  }

  void transposeS() {
    W temp = 0.0;
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
    int info=0;
    /*
     * [  I_{K*K}   A            0       ]  [ x_K ] = [ d_K ]  bandwidth  envery
     *demand
     * [  B         E            0       ]  [ x_N ] = [ c_N ]  capacity of
     *status links
     * [  C         D            I_{J*J} ]  [ x_J ] = [ c_J ]  capacity of other
     *links
     * x_N=( B d_K -c_N )/(BA-E)
     *set S=BA-E, b= B d_K -c_N
     *x_N = b/S
     */
    if (N > 0) {
      fill(b, b + N, 0.0);

      for (int i = 0; i < N; i++) {
        b[i] = -rhs[K + status_links[i]];
      }
      for (int i = 0; i < N; i++) {
        int link = status_links[i];

        const set<int>& pps = status_primary_path_locs[link];
        for (set<int>::const_iterator it = pps.begin(); it != pps.end(); it++) {
          b[i] += rhs[paths[*it].owner];
        }
      }
      copy(S, S + N * N, workS);
      
      if( sizeof( W )==sizeof( double ) ){
        dgesv_(&N, &nrhs, (double*)workS, &lda, ipiv, (double*)b, &ldb, &info);
        
      }else{
        sgesv_(&N, &nrhs, (float*)workS, &lda, ipiv, (float*)b, &ldb, &info);        
      }

      if (info > 0) {
        assert( false );
        printf(
            "The diagonal element of the triangular factor of "
            "A,\n");
        printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
      memcpy(x_N, b, N * sizeof(W));
    }

    /**
     * x_K=d_K-A x_N
     *
     */
    for (int i = 0; i < K; i++) {
      x_K[i] = rhs[i];
    }

    for (int i = 0; i < K; i++) {
      const set<int>& pathindices = demand_second_path_locs[i];
      for (set<int>::const_iterator it = pathindices.begin();
           it != pathindices.end(); it++) {
        x_K[i] -= x_N[getNindex(paths[*it].link) - K];
      }
    }


  }

  EXIT_VARIABLE computeExitVariable() {
    /**
     *  choose enter base as i which a[ i]>0 and x[ i ]/a[ i ]= min x/a
     *
     */

    int reK = -1;
    W minK = numeric_limits<W>::max();
    W temp;
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
    W minN = numeric_limits<W>::max();
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
    W minJ = numeric_limits<W>::max();

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
      if (minN <= EPS) sdata.empty_iterator_num++;
      return re;
    }
    re.id = un_status_links[reJ];
    re.type = OTHER_LINK;
    if (minJ <= EPS) sdata.empty_iterator_num++;
    return re;
  }

 public:
  CG(const G& g, const vector<W>& ws, const vector<C>& caps,
     const vector<Demand>& ds)
      : graph(g), demands(ds), orignal_weights(ws), orignal_caps(caps) ,thread_num( 1 ){

#ifdef _OPENMP
#pragma omp parallel for
    for( int i=0; i<1000; i++ ){
      thread_num=omp_get_num_threads(  );
    }
    
    
#endif // _OPENMP
    info = 0;
    origLink_num = graph.getLink_num();

    CZERO = ((C)1e-6);
    EPS = ((W)1e-6);

    K = demands.size();
    demand_second_path_locs.resize(K);

    A = NULL;

    X = NULL;

    ipiv = NULL;

    b = NULL;

    S = NULL;
    workS = NULL;
    S_maxdim = 0;


    totalB = 0;
    for (typename vector<Demand>::const_iterator it = demands.begin();
         it != demands.end(); it++) {
      totalB += it->bandwidth;
    }

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

  bool is_status_link(const int link) const {
    return find(status_links.begin(), status_links.end(), link) !=
        status_links.end();
  }
  
  void setStatusLink(const int link, const int pid  ){
    paths[ pid ].link=link;
    status_link_path_loc[link]=pid;
  }

  void deletePrimaryPath( const int commodityID ){
    
    int primary_pid=primary_path_loc[ commodityID ];
    const vector<int>& orig_path = paths[primary_pid].path;
    for (vector<int>::const_iterator it = orig_path.begin();
         it != orig_path.end(); it++) {
      if (find(status_links.begin(), status_links.end(), *it) !=
          status_links.end()) {
        status_primary_path_locs[*it].erase(primary_pid);
      }
    }
  }
  void addPrimaryPath(const int commodityID   ){
    int pid=primary_path_loc[commodityID];
    const vector<int>& new_path =
        paths[pid].path;
    for (vector<int>::const_iterator it = new_path.begin();
         it != new_path.end(); it++) {
      if (find(status_links.begin(), status_links.end(), *it) !=
          status_links.end()) {
        status_primary_path_locs[*it].insert(pid);
      }
    }
  }
  
  void addStatusLink( const int link ){
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      if (find(paths[pid].path.begin(), paths[pid].path.end(), link) !=
          paths[pid].path.end()) {
        status_primary_path_locs[link].insert(pid);
      }
    }    
  }
  void deleteStatusLink( const int link ){
      vector<int>::iterator it =
          find(status_links.begin(), status_links.end(), link);
      assert(it != status_links.end());

      status_links.erase(it);
      
      int spid = status_link_path_loc[link];
      paths[spid].link = -1;

      status_link_path_loc.erase(link);

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
    assert(find(un_status_links.begin(), un_status_links.end(), i) !=
           un_status_links.end());
    vector<int>::const_iterator it =
        find(un_status_links.begin(), un_status_links.end(), i);
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
    assert(find(status_links.begin(), status_links.end(), i) !=
           status_links.end());
    vector<int>::const_iterator it =
        find(status_links.begin(), status_links.end(), i);
    return (it - status_links.begin()) + K;
  }

  W getOrigCost(const vector<int>& path) const {
    W re = 0.0;
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
  // bool solveLP(const W * M, const int n, W * bb   ){
  //   if( sizeof( W )==sizeof( double ) ){
  //     vector<int> AAp;
  //     vector<int> AAi;
  //     vector<W> AAx;
  //     int nz=0;
  //     AAp.push_back( 0 );
  //     for( int i=0; i< n; i++ ){
  //       for( int j=0; j< n; j++ ){
  //         if( fabs( M[ i*n+j ])>EPS  ){
  //           nz++;
  //           AAi.push_back( j );
  //           AAx.push_back( M[ i*n+j ] );
  //         }
  //       }
  //       AAp.push_back( nz );
  //     }

  //     klu_symbolic *Symbolic ;
  //     klu_numeric *Numeric ;
  //     klu_common Common ;

  //     int *Ap=&AAp[ 0 ];
  //     int *Ai=&AAi[ 0 ];
  //     double *Ax=&AAx[ 0 ];
      
  //     Symbolic = klu_analyze (n, Ap, Ai, &Common) ;
  //     Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
  //     klu_solve (Symbolic, Numeric, n, 1, b, &Common) ;
  //     klu_free_symbolic (&Symbolic, &Common) ;
  //     klu_free_numeric (&Numeric, &Common) ;
      
  //   }else{
      
  //   }
   
  // }

  bool solve() {
    initial_solution();

    iteration();

    if (info > 0) {
      printResult();
    }
    return true;
  }

  void initial_solution() {
    paths.resize(K);

    primary_path_loc.resize(K, 0);
    vector<C> temp_p(K);

    inf_weight = 0;
    for (int i = 0; i < origLink_num; i++) {
      inf_weight += orignal_weights[i];
    }
    inf_weight *= 4;
    inf_weight += 1;

    vector<int> srcs, snks;
    int src, snk=0;
    for (int i = 0; i < origLink_num; i++) {
      graph.findSrcSnk(i, src, snk);
      srcs.push_back(src);
      snks.push_back(snk);
    }
    
    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      srcs.push_back(src);
      snks.push_back(snk);
      orignal_weights.push_back(inf_weight / 2);
      orignal_caps.push_back(bw);      
    }

    vector<W> temp_cap(orignal_caps);


    vector<bool> success( K, false );
    
    for (int i = 0; i < K; i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      vector<W> ws = orignal_weights;
      for (size_t j = 0; j < origLink_num; j++) {
        if (temp_cap[j] < bw) {
          ws[j] = inf_weight;
        }
      }
      
      vector<int> path;
      if (bidijkstra_shortest_path(graph, ws, inf_weight, src, snk, path)) {
        
        success[ i ]=true;

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
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwidth;
      srcs.push_back(src);
      snks.push_back(snk);
      orignal_weights.push_back(inf_weight / 2);
      orignal_caps.push_back(bw);
      if( !success[ i ]){
        vector<int> path;
        path.push_back( srcs.size(  )-1 );
        paths[i].path = path;
        paths[i].owner = i;
        primary_path_loc[i] = i;
        temp_p[i] = bw;
      }
      success[ i ]=true;
    }


    newGraph.initial(srcs, snks, orignal_weights);

    update_weights = orignal_weights;
    update_caps = orignal_caps;
    N = 0;
    J = srcs.size();
    inf_weight = numeric_limits<W>::max() / 3;

    for (int i = 0; i < J; i++) {
      un_status_links.push_back(i);
    }

    rhs.resize(demands.size() + srcs.size(), (C)0.0);
    for (size_t i = 0; i < K; i++) {
      rhs[i] = demands[i].bandwidth;
    }

    for (size_t i = 0; i < update_caps.size(); i++) {
      rhs[i + K] = update_caps[i];
    }

    dual_solution.resize(J, 0);
    b = new W[J];
    A = new W[K + J];
    X = new W[K + J];
    fill( X, X+J, 0.0 );
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
        if(sdata.iterator_num%10==0  ){
        
          C sobj=success_obj(  );
          std::cout << sdata.iterator_num <<" status link num: "<<N<<  " objvalue: " << computeOBJ()<<" success fractional bw: "<<sobj<<" success rat: "
                    <<sobj/(totalB+0.01  )<<"  the obj gap from opt is: "<<sdata.estimee_opt_diff<< std::endl;
        }
      }


      /**
       *  enter variable choose
       *
       */
      update_edge_left_bandwith();

      ENTER_VARIABLE enter_commodity = chooseEnterPath();
      if (enter_commodity.id < 0) {
        return;
      }

      EXIT_VARIABLE exit_base;

      if (PATH_T == enter_commodity.type) {
        /**
         *  exit base  choose
         *
         */
        exit_base = getExitBasebyPath(enter_commodity);
      } else {
        exit_base = getExitBasebyStatusLink(enter_commodity);
      }

      devote(enter_commodity, exit_base);

      N = status_links.size();
      J = newGraph.getLink_num() - N;
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
   * choose a comodity which the current best path has biggest diff from old
   *solution
   *
   * @return
   */
  ENTER_VARIABLE chooseEnterPath() {
    ENTER_VARIABLE enter_variable;
    enter_variable.type = PATH_T;
    enter_variable.id = -1;
    /**
     *  check status link dual value
     *
     */
    W min_diff = -EPS;
    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      if (dual_solution[link] < min_diff) {
        min_diff = dual_solution[link];
        enter_variable.id = link;
        enter_variable.type = LINK_T;
        enter_variable.path.clear();
      }
    }

    if (enter_variable.id >= 0) {
      return enter_variable;
    }
    sdata.estimee_opt_diff=0;
    vector<W> opt_gap( thread_num, 0 );
    vector<W> min_diffs( thread_num, -EPS );
    vector<vector<int>> path( thread_num );
    vector<ENTER_VARIABLE> enter_variables( thread_num );
    for (int i=0; i < thread_num; i++) {
      enter_variables[ i ].type = PATH_T;
      enter_variables[ i ].id = -1;
    }
#pragma omp parallel for
    for (int i = 0; i < K; i++) {
      int tid= omp_get_thread_num( );

      int src = demands[i].src;
      int snk = demands[i].snk;

      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]].path, (W)0.0);

      if (bidijkstra_shortest_path(newGraph, update_weights, inf_weight, src,
                                   snk, path[ tid ])) {
        W new_cost = path_cost(update_weights, path[ tid ], (W)0.0);
        
        W temp_diff = (new_cost - old_cost);
        if( temp_diff<-EPS ){
          opt_gap[ tid ]-=temp_diff*demands[ i ].bandwidth;

          temp_diff       *=leftBandwith(path[ tid ]  ) +EPS ;
        
          if (temp_diff < min_diffs[ tid ]) {
            min_diffs[ tid ] = temp_diff;

            enter_variables[ tid ].id = i;
            enter_variables[ tid ].path = path[ tid ];
          }          
        }
      }
    }

    int chid=0;
    min_diff=min_diffs[ 0 ];
    for( int i=1; i< thread_num; i++  ){
      if(min_diffs[ i ]<min_diff  ){
        min_diff=min_diffs[ i ];
        chid=i;
      }
    }
    for( int i=0;i< thread_num; i++ ){
      sdata.estimee_opt_diff+=opt_gap[ i ];
    }

    sort(enter_variables[ chid ].path.begin(), enter_variables[ chid ].path.end());

    return enter_variables[ chid ];
    
  }

  /**
   * [  I_{K*K}   A            0       ]  [ a_K ] = [ b_K ]
   * [  B         E            0       ]  [ a_N ] = [ b_N ]
   * [  C         D            I_{J*J} ]  [ a_J ] = [ b_J ]
   * a_N=(B b_K-b_N) /( BA-I )
   *
   * [  I_{K*K}   A            0       ]  [ x_K ] = [ d_K ]  bandwidth  envery
   *demand
   * [  B         I_{N*N}      0       ]  [ x_N ] = [ c_N ]  capacity of status
   *links
   * [  C         D            I_{J*J} ]  [ x_J ] = [ c_J ]  capacity of other
   *links
   * y_N=( B d_K -c_N )/(BA-I)
   * @param k
   *
   * @return
   */

  EXIT_VARIABLE getExitBasebyPath(const ENTER_VARIABLE& enterCommodity) {
    /**
     * A a=d
     *
     */

    fill(A, A + K + N + J, 0.0);
    a_K = A;
    a_N = A + K;
    a_J = A + K + N;

    const vector<int>& path = enterCommodity.path;

    int commodity_primary_path_loc = primary_path_loc[enterCommodity.id];

    const vector<int>& commodity_path = paths[commodity_primary_path_loc].path;

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
       * [  I_{K*K}   A            0       ]  [ a_K ] = [ b_K ]
       * [  B         I_{N*N}      0       ]  [ a_N ] = [ b_N ]
       * [  C         D            I_{J*J} ]  [ a_J ] = [ b_J ]
       * a_N=(B b_K-b_N) /( BA-I )
       */
      int nrhs = 1;
      int lda = N;

      int ldb = N;
      int info=0;

      /**
       * b=B b_K -b_N
       *
       */

      fill(b, b + N, 0.0);
      for (vector<int>::const_iterator it = commodity_path.begin();
           it != commodity_path.end(); it++) {
        vector<int>::iterator fid =
            find(status_links.begin(), status_links.end(), *it);
        if (fid != status_links.end()) {
          b[fid - status_links.begin()] = 1.0;
        }
      }

      for (vector<int>::const_iterator it = path.begin(); it != path.end();
           it++) {
        vector<int>::iterator fid =
            find(status_links.begin(), status_links.end(), *it);
        if (fid != status_links.end()) {
          b[fid - status_links.begin()] -= 1.0;
        }
      }
      /**
       * a_N=( Bb_K-b_N)/( BA-I )=b/S
       *
       */

      copy(S, S + N * N, workS);
      if( sizeof( W )==sizeof( double ) ){
        char c='N';
        // dgetrs_(&c, &N, &nrhs, (double*)workS, &lda, ipiv, (double*)b, &ldb, &info);
        dgesv_(&N, &nrhs, (double*)workS, &lda, ipiv, (double*)b, &ldb, &info);
        
      }else{
        char c='N';
         sgesv_(&N, &nrhs, (float*)workS, &lda, ipiv, (float*)b, &ldb, &info);
        // sgetrs_(&c, &N, &nrhs, (float*)workS, &lda, ipiv, (float*)b, &ldb, &info);        
      }

      if (info > 0) {
        assert( false );
        printf(
            "The diagonal element of the triangular factor of "
            "A,\n");
        printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
      memcpy(a_N, b, N * sizeof(double));

      /**
       * a_K=b_K-A a_N
       *
       */

      a_K[enterCommodity.id] = 1.0;
      for (int i = 0; i < K; i++) {
        for (set<int>::iterator it = demand_second_path_locs[i].begin();
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
        if (fabs(a_K[i])>EPS) {
          left_AK.push_back(i);
        }
      }

      for (size_t i = 0; i < un_status_links.size(); i++) {
        for (int k = 0; k < left_AK.size(); k++) {
          int pid = left_AK[k];
          int ppid = primary_path_loc[pid];
          if (find(paths[ppid].path.begin(), paths[ppid].path.end(),
                   un_status_links[i]) != paths[ppid].path.end())
            a_J[i] -= a_K[pid];
        }
      }

      vector<int> left_AN;

      for (int i = 0; i < N; i++) {
        if (fabs(a_N[i]) >EPS) {
          left_AN.push_back(i);
        }
      }

      for (size_t i = 0; i < un_status_links.size(); i++) {
        for (int k = 0; k < left_AN.size(); k++) {
          int pid = left_AN[k];
          int ppid = status_link_path_loc[status_links[pid]];
          if (find(paths[ppid].path.begin(), paths[ppid].path.end(),
                   un_status_links[i]) != paths[ppid].path.end())
            a_J[i] -= a_N[pid];
        }
      }
    }

    return computeExitVariable();
  }

  EXIT_VARIABLE getExitBasebyStatusLink(const ENTER_VARIABLE& enterLink) {
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
    int info=0;

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
    copy(S, S + N * N, workS);
    if( sizeof( W )==sizeof( double ) ){
      char c='N';
      
      // dgetrs_(&c, &N, &nrhs, (double*)workS, &lda, ipiv, (double*)b, &ldb, &info);        
      dgesv_(&N, &nrhs, (double*)workS, &lda, ipiv, (double*)b, &ldb, &info);        
    }else{
      char c='N';
       sgesv_(&N, &nrhs, (float*)workS, &lda, ipiv, (float*)b, &ldb, &info);
      // sgetrs_(&c, &N, &nrhs, (float*)workS, &lda, ipiv, (float*)b, &ldb, &info);        
    }

    if (info > 0) {
      assert( false );
      printf(
          "The diagonal element of the triangular factor of "
          "A,\n");
      printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
      printf("the solution could not be computed.\n");
      exit(1);
    }
    memcpy(a_N, b, N * sizeof(W));

    /**
     * a_K=-A a_N
     *
     */

    for (int i = 0; i < K; i++) {
      for (set<int>::iterator it = demand_second_path_locs[i].begin();
           it != demand_second_path_locs[i].end(); it++) {
        int pindex = *it;
        int link = paths[pindex].link;
        assert( link>=0 );
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
      for (int k = 0; k < left_AK.size(); k++) {
        int pid = left_AK[k];
        int ppid = primary_path_loc[pid];
        
        if (!paths[ppid].path.empty()&&find(paths[ppid].path.begin(), paths[ppid].path.end(),
                 un_status_links[i]) != paths[ppid].path.end())
          a_J[i] -= a_K[pid];
      }
    }

    vector<int> left_AN;

    for (int i = 0; i < N; i++) {
      if (fabs(a_N[i]) >EPS) {
        left_AN.push_back(i);
      }
    }

    for (size_t i = 0; i < un_status_links.size(); i++) {
      for (int k = 0; k < left_AN.size(); k++) {
        int pid = left_AN[k];
        int ppid = status_link_path_loc[status_links[pid]];
        if (!paths[ppid].path.empty()&&find(paths[ppid].path.begin(), paths[ppid].path.end(),
                 un_status_links[i]) != paths[ppid].path.end())
          a_J[i] -= a_N[pid];
      }
    }

    return computeExitVariable();
  }

  void devote(ENTER_VARIABLE& enter_commodity, EXIT_VARIABLE& exit_base) {
    // std::cout << "devote" << std::endl;

    if (PATH_T == enter_commodity.type) {
      // std::cout << "path" << std::endl;
      if (DEMAND_T == exit_base.type) {
        // std::cout << "demand" << std::endl;
        int exit_commodity_id = exit_base.id;
        int exit_primary_pid = primary_path_loc[exit_commodity_id];

        /**
         * exit primary will been deleted from base matrix
         *
         */
        deletePrimaryPath( exit_commodity_id );


        /**
         * when enter commodity and exit commodity are same then replace the
         *commodity primary path with  enter path
         *
         */

        paths[exit_primary_pid].path = enter_commodity.path;
        paths[exit_primary_pid].owner = enter_commodity.id;

        if (enter_commodity.id != exit_base.id) {
          /**
           * when enter commodity and exit commodity are diff then replace the
           *commodity primary path with the second path of the commodity which
           * crossponding status link cross enter path and make enter path as
           *path coreesponding status link
           *
           */

          assert(  !demand_second_path_locs[exit_commodity_id].empty(  ));
          set<int>::const_iterator it =
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
          paths[pid].link=-1;
          
          demand_second_path_locs[enter_commodity.id].insert(exit_primary_pid);

          setStatusLink(link,  exit_primary_pid );

        }

        addPrimaryPath( exit_commodity_id );

      } else if (STATUS_LINK == exit_base.type) {
        // std::cout << "status link" << std::endl;

        int pid = status_link_path_loc[exit_base.id];
        paths[pid].path = enter_commodity.path;

        /**
         * when the owner of the exit path is not owner of enter path
         *
         */

        if (enter_commodity.id != paths[pid].owner) {
          demand_second_path_locs[paths[pid].owner].erase(pid);
          demand_second_path_locs[enter_commodity.id].insert(pid);
        }
        paths[pid].owner = enter_commodity.id;

      } else {
        
        // std::cout << "other link" << std::endl;
        // exit is un status link then from un_status link to status link
        if (empty_paths.empty()) {
          Path npath(enter_commodity.path, enter_commodity.id, exit_base.id);

          paths.push_back(npath);

          status_links.push_back(exit_base.id);
          status_link_path_loc[exit_base.id] = paths.size() - 1;

          demand_second_path_locs[enter_commodity.id].insert(paths.size() - 1);

        } else {
          int pid = empty_paths.back();
          empty_paths.pop_back();

          paths[pid].path = enter_commodity.path;
          paths[pid].owner = enter_commodity.id;
          paths[pid].link = exit_base.id;

          status_links.push_back(exit_base.id);

          status_link_path_loc[exit_base.id] = pid;
          demand_second_path_locs[enter_commodity.id].insert(pid);
        }

        int link = exit_base.id;
        addStatusLink( link );

      }

    } else {

      // std::cout << "link" << std::endl;
      /**
       * enter a status link
       *
       */
      int enter_status_link = enter_commodity.id;
      
      int spid = status_link_path_loc[enter_status_link];
      deleteStatusLink( enter_status_link );
      

      if (DEMAND_T == exit_base.type) {
        // std::cout << "demand " << std::endl;
        int exit_commodity_id = exit_base.id;
        
        deletePrimaryPath( exit_commodity_id );
        
        if (paths[spid].owner == exit_base.id) {
          primary_path_loc[exit_base.id] = spid;
          demand_second_path_locs[exit_base.id].erase(spid);
        } else {

          assert( !demand_second_path_locs[exit_base.id].empty(  ) );
          set<int>::const_iterator it =
              demand_second_path_locs[exit_base.id].begin();

          int pid = *it;
          int link = paths[pid].link;

          assert(link >= 0);

          demand_second_path_locs[exit_base.id].erase(it);
          primary_path_loc[exit_base.id] = pid;

          setStatusLink( link, spid );
        }
        
        addPrimaryPath( exit_commodity_id );
      } else if (STATUS_LINK == exit_base.type) {

        // std::cout << "status link sssssssssssssssssssssssssssssssssss" << std::endl;
        
        if(exit_base.id == enter_commodity.id  ){
          // std::cout << "status link nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn" << std::endl;
          empty_paths.push_back(spid);
          demand_second_path_locs[paths[spid].owner].erase(spid);
          paths[spid].owner = -1;          
        }else{
          int pid=status_link_path_loc[exit_base.id  ];
          empty_paths.push_back( pid );
          demand_second_path_locs[ paths[ pid ].owner ].erase( pid );
          
          paths[ pid ].owner=-1;
          paths[ pid ].link=-1;
          
          setStatusLink( exit_base.id, spid );
        }

        
      } else {
        // std::cout << "other link " << std::endl;
        status_links.push_back(exit_base.id);
        
        int link = exit_base.id;
        
        setStatusLink( link, spid );
        
        addStatusLink( link );

      }
    }
  }

  void update_edge_left_bandwith() {
    edgeLeftBandwith = update_caps;
    vector<C> allow( K, 0.0 );

    for (int i = 0; i < K; i++) {

      int pid = primary_path_loc[i];
      allow[paths[pid].owner ]+=X[i];
      const vector<int>& path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i];
      }
    }

    for (int i = 0; i < N; i++) {
      int link = status_links[i];
      int pid = status_link_path_loc[link];
      allow[paths[pid].owner ]+=X[i + K];
      const vector<int>& path = paths[pid].path;
      for (vector<int>::const_iterator lit = path.begin(); lit != path.end();
           lit++) {
        edgeLeftBandwith[*lit] -= X[i + K];
      }
    }
    
    for( int i=0; i< N+J; i++ ){
      assert(edgeLeftBandwith[ i ]> -EPS  );
    }
    for( int i=0; i< K; i++ ){
      assert( allow[ i ]-demands[ i ].bandwidth >-EPS && allow[ i ]-demands[ i ].bandwidth <EPS );
    }
    
    x_J=X+N+K;
    for( int i=0; i<J; i++ ){
      int link=un_status_links[ i ];
      x_J[ i ]=edgeLeftBandwith[ link ];
      assert( x_J[ i ]>-EPS );
    }
  }

  C leftBandwith(const vector<int> path) {
    if (path.empty()) return 0;
    C re = edgeLeftBandwith[path.front()];
    for (vector<int>::const_iterator it = path.begin(); it != path.end();
         it++) {
      re = min(re, edgeLeftBandwith[*it]);
    }

    return re;
  }

  /**
   *[ y_K  y_N  0 ]  [  I_{K*K}   A            0       ]   =  [ C_K C_N 0 ]
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
    int info=0;

    fill(b, b + N, 0.0);
    for (int i = 0; i < N; i++) {
      int linkid = status_links[i];
      int pid = status_link_path_loc[linkid];
      b[i] =  -getOrigCost ( paths[ pid ].path );
    }

    fill(A, A + K, 0.0);
    for (int i = 0; i < K; i++) {
      int pid = primary_path_loc[i];
      A[i] = getOrigCost ( paths[ pid ].path );
    }
    for (int i = 0; i < N; i++) {
      int linkid = status_links[i];
      int pid = status_link_path_loc[linkid];
      int oindex = paths[pid].owner;
      b[i] += A[oindex];
    }
    copy(S, S + N * N, workS);
    if( sizeof( W )==sizeof( double ) ){
      dgesv_(&N, &nrhs, (double*)workS, &lda, ipiv, (double*)b, &ldb, &info);        
    }else{
      sgesv_(&N, &nrhs, (float*)workS, &lda, ipiv, (float*)b, &ldb, &info);        
    }

    if (info > 0) {
      assert( false );
      printf(
          "The diagonal element of the triangular factor of "
          "A,\n");
      printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
      printf("the solution could not be computed.\n");
      exit(1);
    }

    update_weights=orignal_weights;
    fill(dual_solution.begin(  ), dual_solution.end(  ), 0.0);

    for (int i = 0; i < N; i++) {
      dual_solution[status_links[i]] =- b[i];
      update_weights[status_links[i]] -= b[i];
    }
  }

  void printResult() {
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
