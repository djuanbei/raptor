
/**
 * @file   mcfgc.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Mon Apr  4 14:36:13 2016
 *
 * @brief  column generation method for multi-commodity flow problem
 *
 *
 */
#include <vector>
#include <limits>
#include <set>
#include <map>

#include <algorithm>

#include "graphalg.hpp"
using namespace fast_graph;
using namespace std;

/* DGESV prototype */
extern "C" {
void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b,
            int* ldb, int* info);
}

namespace mcmcf {

template <typename G, typename W, typename C>
class CG {
 public:
  struct Demand {
    int src, snk;
    C bandwith;
    Demand() : src(0), snk(0), bandwith(0) {}
  };

  struct Flows {
    vector<vector<int> > paths;
    vector<C> bandwiths;
  };

 private:
  G graph;
  int origLink_num;
  G newGraph;
  W inif_weight;

  vector<Demand> demands;

  vector<W> orignal_weights;
  vector<C> orignal_caps;

  vector<W> update_weights;
  vector<C> update_caps;
  vector<C> edgeLeftBandwith;
  vector<vector<int> > paths;
  vector<int> owner;

  vector<C> primal_solution;
  map<int, C> dual_solution;

  vector<int> status_links;

  vector<set<int> > demand_sec_paths;

  vector<set<int> > status_primary_paths;

  vector<int> primary_paths;

  map<int, int> second_paths;

  C CZERO;

 public:
  CG(const G& g, const vector<W>& ws, const vector<C>& caps,
     const vector<Demand>& ds)
      : graph(g), demands(ds), orignal_weights(ws), orignal_caps(caps) {
    origLink_num = graph.getLink_num();

    CZERO = ((C)1e-6);
  }

  int getIndex( int i ) const{
    int re=i;
    for( vector<int>::const_iterator it=status_links.begin(  ); it!= status_links.end(  ); it++ ){
      if( *it>i )  re++;
    }
    return re;
    
  }

  bool solve() {
    if (!initial_solution()) {
      update_edge_left_bandwith();
      iteration();
    }
  }

  bool initial_solution() {
    paths.resize(demands.size());
    primal_solution.resize(demands.size());

    vector<bool> succ_sate(demands.size(), false);
    vector<double> temp_cap(orignal_caps);
    inif_weight = 0;
    for (int i = 0; i < origLink_num; i++) {
      inif_weight += orignal_weights[i];
    }
    inif_weight *= 2;
    inif_weight += 1;

    for (size_t i = 0; i < demands.size(); i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      C bw = demands[i].bandwith;
      vector<W> ws = orignal_weights;
      for (size_t j = 0; j < origLink_num; j++) {
        if (temp_cap[j] < bw) {
          ws[j] = inif_weight;
        }
      }
      vector<int> path;
      if (bidijkstra_shortest_path(graph, ws, inif_weight, src, snk, path)) {
        succ_sate[i] = true;
        for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
          temp_cap[*it] -= bw;
        }
        paths[i] = path;
        primal_solution[i] = bw;
      }
    }
    vector<int> srcs, snks;
    int src, snk;
    for (int i = 0; i < origLink_num; i++) {
      graph.findSrcSnk(i, src, snk);
      srcs.push_back(src);
      snks.push_back(snk);
    }
    update_weights = orignal_weights;
    update_caps = orignal_caps;
    bool re = true;

    for (size_t i = 0; i < demands.size(); i++) {
      if (!succ_sate[i]) {
        re = false;
        int src = demands[i].src;
        int snk = demands[i].snk;
        C bw = demands[i].bandwith;
        srcs.push_back(src);
        snks.push_back(snk);
        update_weights.push_back(inif_weight / 2);
        update_caps.push_back(bw);
        vector<int> path;
        path.push_back(srcs.size() - 1);
        paths[i] = path;
        primal_solution[i] = bw;
      }
    }
    if (!re) {
      newGraph.initial(srcs, snks, update_weights);
      status_links.resize(srcs.size(), -1);  // no status links
    }
    return re;
  }

  void iteration() {
    while (true) {
      int enter_commodity = chooseEnterPath();
      if (enter_commodity < 0) {
        return;
      }

      int exit_base = getExitBase(enter_commodity);

      devote(enter_commodity, exit_base);

      update_edge_left_bandwith();

      update_edge_cost();
    }
  }
  /**
   *
   * choose a comodity which the current best path has biggest diff from old
   *solution
   *
   * @return
   */
  int chooseEnterPath() {
    int enter_commodity = -1;
    W best_diff = 0;
    vector<int> path;
    for (size_t i = 0; i < demands.size(); i++) {
      int src = demands[i].src;
      int snk = demands[i].snk;
      W old_cost = path_cost(update_weights, paths[i], ((W)0.0));

      if (bidijkstra_shortest_path(newGraph, update_weights, inif_weight, src,
                                   snk, path)) {
        W new_cost = path_cost(update_weights, path, ((W)0.0));
        W temp_diff = old_cost - new_cost;
        if (temp_diff > best_diff) {
          best_diff = temp_diff;
          enter_commodity = i;
        }
      }
    }

    return enter_commodity;
  }

  /**
   * [  I_{K*K}   A            0       ]  [ a_K ] = [ b_K ]
   * [  B         I_{N*N}      0       ]  [ a_N ] = [ b_N ]
   * [  C         D            I_{J*J} ]  [ a_J ] = [ b_J ]
   * a_N=(B b_K-b_N) /( BA-I )
   * @param k
   *
   * @return
   */

  int getExitBase(const int k) {
    vector<int> path;
    int src = demands[k].src;
    int snk = demands[k].snk;
    bidijkstra_shortest_path(newGraph, update_weights, inif_weight, src, snk,
                             path);
    int K = demands.size();

    int n = status_links.size();
    if (n > 0) {
      int nrhs = 1;
      int lda = n;
      int* ipiv = new int[n];
      int ldb = n, info;

      /**
       * S=BA-I
       *
       */

      double* S = new double[n * n];

      fill(S, S + n * n, 0.0);
      for (int i = 0; i < n; i++) {
        S[i * n + i] = -1.0;
      }
      for (int i = 0; i < n; i++) {
        if (status_primary_paths[i].empty()) continue;

        for (int j = 0; j < n; j++) {
          int oindex = owner[K + j];
          if (status_primary_paths[i].find(oindex) !=
              status_primary_paths[i].end())
            S[i * n + j] += 1.0;
        }
      }

      /**
       * b=B b_K -b_N
       *
       */
      double* b = new double[n];

      fill(b, b + n, 0.0);
      for (vector<int>::iterator it = paths[k].begin(); it != paths[k].end();
           it++) {
        vector<int>::iterator fid =
            find(status_links.begin(), status_links.end(), *it);
        if (fid != status_links.end()) {
          b[fid - status_links.begin()] = 1.0;
        }
      }

      for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
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

      double* a_N = new double[n];

      dgesv_(&n, &nrhs, S, &lda, ipiv, b, &ldb, &info);
      if (info > 0) {
        printf("The diagonal element of the triangular factor of A,\n");
        printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
        printf("the solution could not be computed.\n");
        exit(1);
      }

      /**
       * a_K=b_K-A a_N
       *
       */

      double* a_K = new double[K];

      fill(a_K, a_K + K, 0.0);
      a_K[k] = 1.0;
      for (int i = 0; i < K; i++) {
        for (set<int>::iterator it = demand_sec_paths[i].begin();
             it != demand_sec_paths[i].end(); it++) {
          a_K[i] -= a_N[*it - K];
        }
      }
      /**
       * a_J=b_J-C a_K-D a_K
       *
       */
      double* a_J=new double [ newGraph.getLink_num(  )-n ];
      fill(a_J,a_J+ newGraph.getLink_num(  )-n , 0.0 );
      for( vector<int>::iterator it=path.begin(  ); it!= path.end(  ); it++  ){
        a_J[ getIndex( *it )-n ]=1.0;
      }

      
      
      delete[] ipiv;
      delete[] S;
      delete[] b;
      delete[] a_N;
      delete[] a_K;
      delete[] a_J;
    }
  }

  void devote(int enter_commodity, int exit_base) {}

  void update_edge_left_bandwith() {
    edgeLeftBandwith = update_caps;
    size_t i = 0;
    for (i = 0; i < paths.size(); i++) {
      for (vector<int>::iterator lit = paths[i].begin(); lit != paths[i].end();
           lit++) {
        edgeLeftBandwith[*lit] -= primal_solution[i];
      }
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

  void update_edge_cost() {}
};
}
