
/**
 * @file   mcfgc.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Mon Apr  4 14:36:13 2016
 *
 * @brief  column generation method for multi-commodity flow problem
 *
 *
 */
#include<deque>
#include <vector>
#include <limits>
#include <set>
#include <map>
#include <cstring>
#include <cassert>
#include <limits>
#include <algorithm>

#include "graphalg.hpp"
using namespace fast_graph;
using namespace std;

/* DGESV prototype */
extern "C" {
void dgesv_(int* N, int* nrhs, double* a, int* lda, int* ipiv,
            double* b, int* ldb, int* info);
}

namespace mcmcf {
enum ENTER_BASE_TYPE{
  PATH_T=0,
  LINK_T=1
  
};

struct ENTER_VARIABLE
{
  ENTER_BASE_TYPE type;
  int id;
  vector<int> path;
  ENTER_VARIABLE(  ):type( PATH_T ),id( 0 ){
    
  }
};

enum EXIT_BASE_TYPE{
  DEMAND_T=0,
  STATUS_LINK=1,
  OTHER_LINK=2
  
};
struct EXIT_VARIABLE
{
  EXIT_BASE_TYPE type;
  int id;
  EXIT_VARIABLE(  ):type( DEMAND_T ), id( 0 ){
    
  }
  
};

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
  
  vector<vector<int> > paths;//all the paths save in this vector
  vector<int> empty_paths; // the location which  is delete path
  vector<int> owner; // every path has over demand
  vector<int> link_of_path;
    
  vector<C> dual_solution;

  vector<bool> is_status_link;
  
  vector<int> status_links; // contain the index of status links
  
  vector<int> un_status_links;

  vector<set<int> > demand_second_path_locs; // belong to fixed demands' second paths

  map<int, set<int> > status_primary_path_locs; // primary paths which corross the status link

  vector<int> primary_path_loc; // every demands have a primary path

  map<int, int> status_link_path_loc; //  the path corresponding status link


  
  vector<C> rhs;
  double *A;

  int * ipiv;
  double *b;

  double *Y;

  double *S;
  double* workS;
  int S_maxdim;

  double* a_K=A;
  double* a_N=A+K;
  double* a_J=A+K+N;
  
  double* y_K;
  double* y_N;
  double* y_J;

  
  C CZERO;
  int K, N, J;
  W EPS;

  void  allocateS( const int N ){
    if(S_maxdim>=N  ) return;
    if( NULL!=S ){
      delete[  ] S;
      S=NULL;
      delete[  ] workS;
    }
    S_maxdim=1.3*N+1;
    S=new double[S_maxdim*S_maxdim  ];
    workS=new double[S_maxdim*S_maxdim  ];
    
  }
  /** 
   * row primary
   * 
   */
  void computeS(  ){
    if( 0==N ) return;
    /**
     * S=BA-I
     *
     */
    if( N>S_maxdim )allocateS( N );
    fill(S, S + N * N, 0.0);
      
    for (int i = 0; i < N; i++) {
      S[i * N + i] = -1.0;
    }
    for (int i = 0; i < N; i++) {
      if (status_primary_path_locs[i].empty()) continue;

      for (int j = 0; j < N; j++) {

        int pindex=status_link_path_loc[ status_links[ j ] ];
        int oindex = owner[pindex];
          
        if (status_primary_path_locs[i].find(oindex) !=
            status_primary_path_locs[i].end())
          S[i * N + j] += 1.0;
      }
    }
  }

  void transposeS(  ){
    double temp=0.0;
    for( int i=0; i< N-1; i++ ){
      for( int j=i+1; j< N; j++ ){
        temp=S[ i*N+j ];
        S[ i*N+j ]= S[ j*N+i ];
        S[ j*N+i ]=temp;
      }
    }
  }

  void  computeRHS(  ){
    int nrhs = 1;
    int lda = N;

    int ldb = N;
    int info;
    /*
     * [  I_{K*K}   A            0       ]  [ y_K ] = [ d_K ]  bandwidth  envery  demand   
     * [  B         I_{N*N}      0       ]  [ y_N ] = [ c_N ]  capacity of status links
     * [  C         D            I_{J*J} ]  [ y_J ] = [ c_J ]  capacity of other links
     * y_N=( B d_K -c_N )/(BA-I)
     *set S=BA-I, b= B d_K -c_N
     *y_N = b/S
     */
    fill( b, b+N, 0.0 );

    for( int i=0; i< N; i++ ){
      b[ i ]=-rhs[ status_links[ i ] ];
    }
    for( int i=0; i< N; i++ ){
      if( status_primary_path_locs.find( status_links[ i ] )!= status_primary_path_locs.end(  )){
        const set<int>& pathindices=status_primary_path_locs.at( i );
        for( set<int>::const_iterator it=pathindices.begin(  ); it!= pathindices.end(  ); it++ ){
          b[ i ]+=rhs[ *it ];
        }
      }
    }
    copy( S, S+N*N, workS );
    dgesv_(&N, &nrhs, workS, &lda, ipiv, b, &ldb, &info);
    if (info > 0) {
      printf(
          "The diagonal element of the triangular factor of "
          "A,\n");
      printf("U(%i,%i) is zero, so that A is singular;\n", info,
             info);
      printf("the solution could not be computed.\n");
      exit(1);
    }
    memcpy( y_N, b, N*sizeof( double )  );



    /**
     * y_K=d_K-A y_N
     *
     */
    for( int i=0; i< K; i++ ){
      y_K[ i ]= rhs[ i ];
    }
      
    for( int i=0; i< K; i++ ){
      const set<int> & pathindices=demand_second_path_locs[ i ];
      for( set<int>::const_iterator it=pathindices.begin(  ); it!= pathindices.end(  ); it++ ){
        y_K[ i ]-=y_N[getNindex(link_of_path[ *it ]  ) -K ];
      }
    }


    

    /**
     * y_J=c_J -C y_K -D y_N
     * 
     */

    for( int i=0 ; i< J; i++ ){
      int linkid=getJIndex( un_status_links[ i ] )-N-K;
      y_J[ linkid ]=rhs[ un_status_links[ i ] ];
    }

    vector<int> left_YK;
    for (int i = 0; i < K; i++) {
      if (y_K[i] != 0.0) {
        left_YK.push_back(i);
      }
    }

    for( size_t i=0; i< un_status_links.size(  ); i++ ){
      for (int k = 0; k < left_YK.size(); k++) {
        int pid = left_YK[k];
        int ppid=primary_path_loc[ pid ];
        if (find(paths[ ppid].begin(), paths[ppid].end(), un_status_links[i]) !=
            paths[ppid].end())
          y_J[i] -= y_K[pid];
      }        
    }

    vector<int> left_YN;

    for (int i = 0; i < N; i++) {
      if (y_N[i] != 0.0) {
        left_YN.push_back(i);
      }
    }

    for( size_t i=0; i< un_status_links.size(  ); i++ ){

      for (int k = 0; k < left_YN.size(); k++) {
        int pid = left_YN[k];
        int ppid= status_link_path_loc[ status_links[ pid ] ];
        if (find(paths[ppid].begin(), paths[ppid].end(), un_status_links[i]) !=
            paths[ppid].end())
          y_J[i] -= y_N[pid];
      }
    }    
    
  }

  EXIT_VARIABLE computeExitVariable(  ){
        /**
     *  choose enter base as i which a[ i]>0 and y[ i ]/a[ i ]= min y/a
     * 
     */


    int reK=-1;
    double minK=numeric_limits<double>::max(  );
    double temp;
      int i=0;
      while (i<K ) {
        if(a_K[ i ]>0  ) {
          reK=i;
          minK=y_K[ i ]/a_K[ i ];
          break;
        }
        i++;
      }
      while(i<K){
        if( a_K[ i ]>0 ){
          temp=y_K[ i ]/a_K[ i ];
          if( temp<minK ){
            reK=i;
            minK=temp;
          }
        }
        i++;
      }

      int reN=-1;
      double minN=numeric_limits<double>::max(  );
      i=0;
      while (i<N ) {
        if(a_N[ i ]>0  ) {
          reN=i;
          minN= y_N[ i ]/a_N[ i ];
          break;
        }
        i++;
      }

      while(i<N){
        if( a_N[ i ]>0 ){
          temp=y_N[i ]/a_N[ i ];
          if( temp<minN ){
            reN=i;
            minN=temp;
          }
        }
        i++;
      }

      

      int reJ=-1;
      double minJ=numeric_limits<double>::max(  );

      i=0;
      while (i<J ) {
        if(a_J[ i ]>0  ) {
          reJ=i;
          minJ=y_J[ i ]/a_J[ i ];
          break;
        }
        i++;
      }

      while(i<J){
        if( a_J[ i ]>0 ){
          temp=y_J[  i  ]/a_J[ i ];
          if( temp<minJ ){
            reJ=i;
            minJ=temp;
          }
        }
        i++;
      }

      EXIT_VARIABLE re;

      if( minK<=minN && minK<=minJ ){
        re.id=reK;
        re.type=DEMAND_T;
        return re;

      }
      if( minN<=minJ ){
        re.id=reN;
        re.type=STATUS_LINK;
        return re;
      }
      re.id=reJ;
      re.type=OTHER_LINK;
      return re;
  }

 public:
  CG(const G& g, const vector<W>& ws, const vector<C>& caps,
     const vector<Demand>& ds)
      : graph(g), demands(ds), orignal_weights(ws), orignal_caps(caps) {
    origLink_num = graph.getLink_num();

    CZERO = ((C)1e-6);
    EPS=( ( W )1e-6 );
    
    K=demands.size(  );

    A=NULL;

    ipiv= new int[ K ];
    b=NULL;
    
    S=NULL;
    workS=NULL;
    S_maxdim=0;
    
    
  }
  ~CG(  ){
    if( NULL!=A ){
      delete[  ] A;
      A=NULL;
          
    }

    delete[  ]  ipiv;
    ipiv=NULL;
    if( NULL!=b ){
      delete[  ] b;
      b=NULL;
    }

    if( NULL!= Y ){
      delete[  ] Y;
      Y=NULL;
    }
    if( NULL !=S ){
      delete[  ] S;
      S=NULL;
      delete[  ] workS;
      workS=NULL;
    }

  }
  /** 
   * 
   * 
   * @param i the id of un_statrus link
   * 
   * @return  the index of the un_status link in simplex matrix
   */
  int getJIndex(int i) const {
    assert(find(status_links.begin(  ), status_links.end(  ), i)==status_links.end(  ) );
    int re = i;
    for (vector<int>::const_iterator it = status_links.begin();
         it != status_links.end(); it++) {
      if (*it > i) re++;
    }
    return re+K;
  }
  /** 
   * 
   * 
   * @param i the id of status link
   * 
   * @return the index of the status link in smplex matrix
   */
  int getNindex( int i ) const{
    assert(find(status_links.begin(  ), status_links.end(  ), i)==status_links.end(  ) );
    int re=0;
    for (vector<int>::const_iterator it = status_links.begin();
         it != status_links.end(); it++) {
      if (*it < i) re++;
    }
    return re+K;
    
  }

  double getOrigCost ( const vector<int> &path ) const{
    double re=0.0;
    for( vector<int>::const_iterator it=path.begin(  ); it!= path.end(  ); it++ ){
      re+=orignal_weights[ *it ];
    }
    return re;
  }

  void computIndexofLinks(  ){
    
    stable_sort(status_links.begin(  ), status_links.end(  )  );
    un_status_links.clear(  );

    size_t    i=0;
    for(int j=0; j< N+J; j++ ){
      if(i>=status_links.size(  )){
        un_status_links.push_back( j );
      }else {
       if (j< status_links[ i ]  ){
         un_status_links.push_back( j );
       }else{
         i++;
       }
      }
    }
  }

  bool solve() {
    initial_solution();

    iteration();

  }

  void initial_solution() {
    
    paths.resize(K);
    owner.resize( K,-1 );
    
    primary_path_loc.resize(K, 0  );
    
    
    vector<bool> succ_sate(K, false);
    vector<double> temp_cap(orignal_caps);
    inif_weight = 0;
    for (int i = 0; i < origLink_num; i++) {
      inif_weight += orignal_weights[i];
    }
    inif_weight *= 2;
    inif_weight += 1;

    for (int i = 0; i < K; i++) {
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
        sort( path.begin(  ), path.end(  ) );
        paths[i] = path;
        
        owner[ i ]=i;
        primary_path_loc[ i ]=i;
        Y[ i]= bw;


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


    for (int i = 0; i < K; i++) {
      if (!succ_sate[i]) {

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
        owner[ i ]=i;
        primary_path_loc[ i ]=i;
        Y[ i]=bw;

        
      }
    }

    N=0;
    J=srcs.size(  );

    is_status_link.resize( J, false );
    for( int i=0; i< J; i++ ){
      un_status_links.push_back( i);
    }

    
    rhs.resize(demands.size()+srcs.size(  ), (C)0.0  );
    for( size_t i=0; i< K; i++ ) {
      rhs[ i ]=demands[ i ].bandwith;
    }
      
    for( size_t i=0; i<update_caps.size(  ); i++  ){
      rhs[ i+K ]=update_caps[ i ];
    }
      
    newGraph.initial(srcs, snks, update_weights);

    // status_link_path_loc.resize( J,-1 );
    
    dual_solution.resize( K+J,0 ) ;
    b=new double[ J ];
    A=new double[ K+J ];
    Y=new double[ K+J ];


  }

  void iteration() {
    while (true) {

      fill( A, A+N+J, 0.0 );
      a_K=A;
      a_N=A+K;
      a_J=A+K+N;
      
      fill( Y, Y+N+J, 0.0 );
      y_K=Y;
      y_N=Y+K;
      y_J=Y+K+N;

      
      update_edge_left_bandwith();
      /**
       *  enter variable choose
       * 
       */

      ENTER_VARIABLE enter_commodity = chooseEnterPath(); 
      if (enter_commodity.id < 0) {
        return;
      }
      EXIT_VARIABLE exit_base;

      if( PATH_T ==enter_commodity.type ){
        /**
         *  exit base  choose
         * 
         */
        exit_base = getExitBasebyPath(enter_commodity);        
      }
      else{
        exit_base = getExitBasebyStatusLink(enter_commodity);
      }
      
      devote(enter_commodity, exit_base);

      N = status_links.size();
      J=newGraph.getLink_num(  )-N;
      computIndexofLinks(  );
      
      computeS(  );
      update_edge_cost();
      
      /**
       * column primary
       * 
       */

      transposeS(  );
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

    W min_diff =0.0;
    vector<int> path;
    ENTER_VARIABLE  enter_variable;
    enter_variable.type=PATH_T;
    enter_variable.id=-1;
    
    for (size_t i = 0; i < demands.size(); i++) {
      
      int src = demands[i].src;
      int snk = demands[i].snk;
      
      W old_cost =
          path_cost(update_weights, paths[primary_path_loc[i]], ((W)0.0));


      if (bidijkstra_shortest_path(newGraph, update_weights, inif_weight, src,
                                   snk, path)) {
        W new_cost = path_cost(update_weights, path, ((W)0.0));
        W temp_diff =  new_cost-old_cost;
        if (temp_diff < min_diff) {
          min_diff = temp_diff;
          
          enter_variable.id = i;
          enter_variable.path=path;
        }
      }
    }
    
    /**
     *  check status link dual value
     * 
     */

    for(int i=0; i< N; i++  ){
      if(dual_solution[ K+i ]< min_diff  ){
        min_diff = dual_solution[ K+i ];
        enter_variable.id = i;
        enter_variable.type=LINK_T;
        enter_variable.path.clear(  );
      }
    }
    sort(enter_variable.path.begin(  ), enter_variable.path.end(  )  );
    
    return enter_variable;
  }

  /**
   * [  I_{K*K}   A            0       ]  [ a_K ] = [ b_K ]
   * [  B         I_{N*N}      0       ]  [ a_N ] = [ b_N ]
   * [  C         D            I_{J*J} ]  [ a_J ] = [ b_J ]
   * a_N=(B b_K-b_N) /( BA-I )
   *
   * [  I_{K*K}   A            0       ]  [ y_K ] = [ d_K ]  bandwidth  envery  demand   
   * [  B         I_{N*N}      0       ]  [ y_N ] = [ c_N ]  capacity of status links
   * [  C         D            I_{J*J} ]  [ y_J ] = [ c_J ]  capacity of other links
   * y_N=( B d_K -c_N )/(BA-I)
   * @param k
   *
   * @return
   */

  EXIT_VARIABLE  getExitBasebyPath(const ENTER_VARIABLE& enterCommodity) {
    
    const vector<int> &path=enterCommodity.path;



    int commodity_primary_path_loc=primary_path_loc[ enterCommodity.id ];

    const vector<int> &commodity_path = paths[ commodity_primary_path_loc ];

    /**
     *  status links are  empty
     * 
     */
    if ( 0==N ){
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

      
      for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
        a_J[getJIndex(*it)-N-K ] = 1.0;
      }


      for (vector<int>::const_iterator it = commodity_path.begin(); it != commodity_path.end(); it++) {
        a_J[getJIndex(*it)-N-K ] -= 1.0;
      }

      /**
       * [ I_{K*K}   0       ] [ y_K ]  =[ rhs_K ]
       * [ B         I_{J*J }] [ y_J ] = [ rhs_J ]
       *
       *  y_K = rhs_K
       *  B rhs_K + y_J = b_J
       *  y_J = b_J -B y_K
       */

      for( int i=0; i< K; i++ ){
        y_K[ i ]=rhs[ i ];
      }
      
      /**
       * y_J = b_J -B y_K
       * 
       */

      for( int i=0; i< J; i++  ){
        y_J[ i ]=edgeLeftBandwith[ un_status_links[ i ] ];
      }

    }
    else {
      
      int nrhs = 1;
      int lda = N;

      int ldb = N;
      int info;

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

      for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
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

      copy( S, S+N*N, workS );
      dgesv_(&N, &nrhs, workS, &lda, ipiv, b, &ldb, &info);
      if (info > 0) {
        printf(
            "The diagonal element of the triangular factor of "
            "A,\n");
        printf("U(%i,%i) is zero, so that A is singular;\n", info,
               info);
        printf("the solution could not be computed.\n");
        exit(1);
      }
      memcpy( a_N, b, N*sizeof( double )  );

      /**
       * a_K=b_K-A a_N
       *
       */

      a_K[enterCommodity.id] = 1.0;
      for (int i = 0; i < K; i++) {
        for (set<int>::iterator it = demand_second_path_locs[i].begin();
             it != demand_second_path_locs[i].end(); it++) {
          int pindex=*it;
          int link=link_of_path[ pindex ];
          a_K[i] -= a_N[getNindex( link ) - K];
        }
      }
      /**
       * a_J=b_J-C a_K-D a_N
       *
       */
      
      for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
        if(is_status_link[ *it ])
          a_J[getJIndex(*it)-N-K ] = 1.0;
      }
      
      vector<int> left_AK;
      for (int i = 0; i < K; i++) {
        if (a_K[i] != 0.0) {
          left_AK.push_back(i);
        }
      }

      for( size_t i=0; i< un_status_links.size(  ); i++ ){
        for (int k = 0; k < left_AK.size(); k++) {
          int pid = left_AK[k];
          int ppid=primary_path_loc[ pid ];
          if (find(paths[ ppid].begin(), paths[ppid].end(), un_status_links[i]) !=
              paths[ppid].end())
            a_J[i] -= a_K[pid];
        }        
      }

      vector<int> left_AN;

      for (int i = 0; i < N; i++) {
        if (a_N[i] != 0.0) {
          left_AN.push_back(i);
        }
      }

      for( size_t i=0; i< un_status_links.size(  ); i++ ){

        for (int k = 0; k < left_AN.size(); k++) {
          int pid = left_AN[k];
          int ppid= status_link_path_loc[ status_links[ pid ] ];
          if (find(paths[ppid].begin(), paths[ppid].end(), un_status_links[i]) !=
              paths[ppid].end())
            a_J[i] -= a_N[pid];
        }
      }
      
      computeRHS(  );
    }
    

    return computeExitVariable(  );

      
  }

  EXIT_VARIABLE  getExitBasebyStatusLink(const ENTER_VARIABLE& enterLink) {

    int nrhs = 1;
    int lda = N;

    int ldb = N;
    int info;
      
    // computeS(  );
    fill(b, b+N, 0.0  );

    b[getNindex( enterLink.id )-K]=-1.0;

    /**
     * a_N=( Bb_K-b_N)/( BA-I )=b/S
     *
     */
    copy( S, S+N*N, workS );
    dgesv_(&N, &nrhs, workS, &lda, ipiv, b, &ldb, &info);
    if (info > 0) {
      printf(
          "The diagonal element of the triangular factor of "
          "A,\n");
      printf("U(%i,%i) is zero, so that A is singular;\n", info,
             info);
      printf("the solution could not be computed.\n");
      exit(1);
    }
    memcpy( a_N, b, N*sizeof( double )  );


    /**
     * a_K=-A a_N
     *
     */

    
    for (int i = 0; i < K; i++) {
      for (set<int>::iterator it = demand_second_path_locs[i].begin();
           it != demand_second_path_locs[i].end(); it++) {
        int pindex=*it;
        int link=link_of_path[ pindex ];
        a_K[i] -= a_N[getNindex( link ) - K];
      }
    }

    /**
     * a_J=-C a_K-D a_N
     *
     */

    vector<int> left_AK;
    for (int i = 0; i < K; i++) {
      if (a_K[i] != 0.0) {
        left_AK.push_back(i);
      }
    }

    for( size_t i=0; i< un_status_links.size(  ); i++ ){
      for (int k = 0; k < left_AK.size(); k++) {
        int pid = left_AK[k];
        int ppid=primary_path_loc[ pid ];
        if (find(paths[ ppid].begin(), paths[ppid].end(), un_status_links[i]) !=
            paths[ppid].end())
          a_J[i] -= a_K[pid];
      }        
    }


    vector<int> left_AN;

    for (int i = 0; i < N; i++) {
      if (a_N[i] != 0.0) {
        left_AN.push_back(i);
      }
    }

    for( size_t i=0; i< un_status_links.size(  ); i++ ){

      for (int k = 0; k < left_AN.size(); k++) {
        int pid = left_AN[k];
        int ppid= status_link_path_loc[ status_links[ pid ] ];
        if (find(paths[ppid].begin(), paths[ppid].end(), un_status_links[i]) !=
            paths[ppid].end())
          a_J[i] -= a_N[pid];
      }
    }

    return computeExitVariable(  ) ;
    
  }

  void devote(ENTER_VARIABLE& enter_commodity,  EXIT_VARIABLE& exit_base) {
    
    if(PATH_T==enter_commodity.type  ){
      
      if(DEMAND_T== exit_base.type ){

        const vector<int>& orig_path=paths[ primary_path_loc[enter_commodity.id  ] ];
        for( vector<int>::const_iterator it=orig_path.begin(  ); it!= orig_path.end(  ); it++ ){
          if(find( status_links.begin(  ), status_links.end(  ), *it )!= status_links.end(  )){
            status_primary_path_locs[ *it ].erase( primary_path_loc[enter_commodity.id  ] );
          }
        }

        /**
         * when enter commodity and exit commodity are same then replace the commodity primary path with  enter path
         * 
         */
        
        paths[ primary_path_loc[enter_commodity.id  ] ]=enter_commodity.path;


        if(enter_commodity.id!=exit_base.id  ){

          /**
           * when enter commodity and exit commodity are diff then replace the commodity primary path with the second path of the commodity which
           * crossponding status link cross enter path and make enter path as path coreesponding status link
           * 
           */
          bool state=false;

          
          for( set<int>::const_iterator it= demand_second_path_locs[ exit_base.id ].begin(  );  it!= demand_second_path_locs[ exit_base.id ].end(  ); it++ ){
            int pid=*it;
            int link=link_of_path[ pid ];

            assert( link>=0 );
            
            if( find(enter_commodity.path.begin(  ), enter_commodity.path.end(  ), link  )!=enter_commodity.path.end(  ) ){
              demand_second_path_locs[ exit_base.id ].erase( it );
              state=true;
              link_of_path[ primary_path_loc[enter_commodity.id  ] ]= link ;
              primary_path_loc[enter_commodity.id  ]=pid;
            }
          }
          assert( state );
        }

        
        const vector<int>& new_path=paths[ primary_path_loc[enter_commodity.id  ] ];
        for( vector<int>::const_iterator it=new_path.begin(  ); it!= new_path.end(  ); it++ ){
          if(find( status_links.begin(  ), status_links.end(  ), *it )!= status_links.end(  )){
            status_primary_path_locs[ *it ].insert( primary_path_loc[enter_commodity.id  ] );
          }
        }
        
      }else if(STATUS_LINK==exit_base.type  ){

        /**
         * the exit status link must cross the enter path
         * 
         */

        assert( find(enter_commodity.path.begin(  ), enter_commodity.path.end(  ),  exit_base.id  )!= enter_commodity.path.end(  ) );
        
        int pid=status_link_path_loc[ exit_base.id ];
        paths[ pid ]=enter_commodity.path;

        /**
         * when the owner of the exit path is not owner of enter path
         * 
         */
        
        if(enter_commodity.id!=  owner[ pid ]  ){
          demand_second_path_locs[owner[ pid ]  ].erase( pid );
          demand_second_path_locs[ enter_commodity.id ].insert( pid );
        }
        owner[ pid ]=enter_commodity.id;
        
      }else {
        // exit is un status link
        if(empty_paths.empty(  )  ){
          paths.push_back( enter_commodity.path );
          owner.push_back( enter_commodity.id );
          link_of_path.push_back(exit_base.id  );
          status_links.push_back( exit_base.id );
          status_link_path_loc[exit_base.id  ]=paths.size(  )-1;
                    
        }else{
          int pid=empty_paths.back(  );
          paths[ pid ]= enter_commodity.path ;
          owner[pid]= enter_commodity.id ;
          link_of_path[ pid ]=exit_base.id;
          empty_paths.pop_back(  );
          status_links.push_back( exit_base.id );
          status_link_path_loc[exit_base.id  ]=pid;
        }
        for( int i=0; i< K; i++ ) {
          int pid=primary_path_loc[ i ];
          if( find(paths[ pid ].begin(  ), paths[ pid ].end(  ),  exit_base.id )!= paths[ pid ].end(  )){
            status_primary_path_locs[ exit_base.id ].insert( pid );
          }
        }
      }
      
    }else{
      /**
       * enter a status link
       * 
       */
      int spid=status_link_path_loc[ enter_commodity.id ];
      status_link_path_loc.erase( enter_commodity.id );
      
      status_primary_path_locs.erase( enter_commodity.id );
      vector<int>::iterator it= find(status_links.begin(), status_links.end(), enter_commodity.id  );
      assert( it!= status_links.end() );
      status_links.erase( it );
      
      is_status_link[enter_commodity.id  ]=false;
      
      
      if(DEMAND_T== exit_base.type ){

        if( owner[ spid ]== exit_base.id ){
          primary_path_loc[exit_base.id  ]=spid;
          demand_second_path_locs[exit_base.id   ].erase( spid );
        }else{

          bool state=false;
          
          
          for( set<int>::const_iterator it= demand_second_path_locs[ exit_base.id ].begin(  );  it!= demand_second_path_locs[ exit_base.id ].end(  ); it++ ){
            int pid=*it;
            int link=link_of_path[ pid ];

            assert( link>=0 );
            
            if( find(paths[ spid ].begin(  ), paths[ spid ].end(  ), link  )!=paths[ pid ].end(  ) ){
              demand_second_path_locs[ exit_base.id ].erase( it );
              state=true;
              link_of_path[ spid ]= link ;
              primary_path_loc[exit_base.id  ]=pid;
            }
          }
          assert( state );
          
        }
        
      }else if(STATUS_LINK==exit_base.type  ){

        assert(exit_base.id== enter_commodity.id  );
                
        
      }else {
        status_links.push_back( exit_base.id );
        is_status_link[exit_base.id  ]=true;
        status_link_path_loc[ exit_base.id ]=spid;

        
      }
      
    }
    
  }
  
  void update_edge_left_bandwith() {
    edgeLeftBandwith = update_caps;

    for( int i=0; i< K; i++ ){
      int pid=primary_path_loc[ i ];
      const vector<int> &path=paths[ pid ];
      for (vector<int>::iterator lit = paths[i].begin(); lit != paths[i].end();
           lit++) {
        edgeLeftBandwith[*lit] -= Y[i];
      }
    }

    for( int i=0; i< N; i++ ){
      int link=status_links[ i ];
      int pid=status_link_path_loc[ link ];
      const vector<int> &path=paths[ pid ];
      for (vector<int>::iterator lit = paths[i].begin(); lit != paths[i].end();
           lit++) {
        edgeLeftBandwith[*lit] -= Y[i];
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

  /**
   *[ Z_K  Z_N  0 ]  [  I_{K*K}   A            0       ]   =  [ C_K C_N 0 ]
   *                 [  B         I_{N*N}      0       ]   
   *                 [  C         D            I_{J*J} ]
   *
   * Z_K+Z_N B=C_K
   * Z_K A + Z_N =C_N
   * Z_N( BA-I ) = C_K A - C_N
   *
   * Z_K= C_K -Z_N B
   */
  void update_edge_cost() {
    int nrhs = 1;
    int lda = N;

    int ldb = N;
    int info;
    
    fill( b, b+N, 0.0 );
    for( int i=0; i< N; i++ ){
      int linkid=status_links[ i ];
      int pid=status_link_path_loc[ linkid ];
      b[ i ]=-getOrigCost ( paths[ pid ] );
    }

    fill( A, A+K, 0.0 ) ;
    for( int i=0; i< K; i++ ){
      int pid=primary_path_loc[ i ];
      A[ i ]=getOrigCost ( paths[ pid ] );
    }
    for( int i=0; i< N; i++ ){
      int linkid=status_links[ i ];
      int pid=status_link_path_loc[ linkid ];
      int oindex=owner[ pid ];
      b[ i]+=A[ oindex ];
    }
    copy( S, S+N*N, workS );
    dgesv_(&N, &nrhs, workS, &lda, ipiv, b, &ldb, &info);
    if (info > 0) {
      printf(
          "The diagonal element of the triangular factor of "
          "A,\n");
      printf("U(%i,%i) is zero, so that A is singular;\n", info,
             info);
      printf("the solution could not be computed.\n");
      exit(1);
    }

    update_weights=orignal_weights;
    fill(dual_solution.begin(  ), dual_solution.end(  ), 0.0);

    for( int i=0; i< N; i++ ){
      dual_solution[ status_links[ i ] ]=-b[ i ];
      update_weights[ status_links[ i ] ]-=b[ i ];
    }
    
  }
};
}
