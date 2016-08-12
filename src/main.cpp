#include "graph.h"
#include "System.h"
#include "graphalg.hpp"
#include <set>
#include <vector>
#include <utility>
#include <cstdlib>
#include <iostream>
#include <random>

#include "mcfcg.hpp"


using namespace std;
using namespace mcmcf;

using namespace fast_graph;

void MCFexample1(  ){
  vector<int> srcs;
  vector<int> snks;
  vector<double> weights, caps;
  srcs.push_back( 0 );
  snks.push_back( 1 );
  weights.push_back( 2 );
  caps.push_back( 1 );

  srcs.push_back( 0 );
  snks.push_back( 1 );
  weights.push_back( 1 );
  caps.push_back( 2 );


  srcs.push_back( 2 );
  snks.push_back( 0 );
  weights.push_back( 1 );
  caps.push_back( 2 );


  srcs.push_back( 2 );
  snks.push_back( 1 );
  weights.push_back( 4 );
  caps.push_back( 2 );


  // srcs.push_back( 0 );
  // snks.push_back( 1 );
  // weights.push_back( 1 );
  // caps.push_back( 2 );
  
  compressed_sparse_row_graph<int,int> graph;

  typedef CG<compressed_sparse_row_graph<int, int>,  double> CG_T;
  vector<CG_T::Demand> demands;
  CG_T::Demand d1;
  d1.src=0;
  d1.snk=1;
  d1.bandwidth=1;
  demands.push_back( d1 );
  CG_T::Demand d2;
  d2.src=2;
  d2.snk=1;
  d2.bandwidth=2;

  demands.push_back( d2 );

  vector<int> ws( weights.size(  ), 1 );

  graph.initial( srcs, snks, ws);
  CG_T cg( graph,weights, caps, demands );
  cg.setInfo( 1 );
  cg.solve(  );
  
      
}

void MCFexample2(  ){

  vector<int> srcs;
  vector<int> snks;
  vector<double> weights, caps;
  srcs.push_back( 0 );
  snks.push_back( 2 );
  weights.push_back( 1 );
  caps.push_back( 1 );

  srcs.push_back( 2 );
  snks.push_back( 1 );
  weights.push_back( 1 );
  caps.push_back( 1 );


  srcs.push_back( 1 );
  snks.push_back( 0 );
  weights.push_back( 1 );
  caps.push_back( 1 );


  
  compressed_sparse_row_graph<int,int> graph;

  typedef CG<compressed_sparse_row_graph<int, int>, double> CG_T;
  vector<CG_T::Demand> demands;
  CG_T::Demand d1;
  d1.src=0;
  d1.snk=1;
  d1.bandwidth=1;
  demands.push_back( d1 );
  
  CG_T::Demand d2;
  d2.src=1;
  d2.snk=2;
  d2.bandwidth=1;

  demands.push_back( d2 );

  CG_T::Demand d3;
  d3.src=2;
  d3.snk=0;
  d3.bandwidth=1;

  demands.push_back( d3 );
  vector<int> ws( weights.size(  ), 1 );
  graph.initial( srcs, snks, ws);
  CG_T cg( graph,weights, caps, demands );
  cg.setInfo( 1 );
  cg.solve(  );
}

void randMCF( const int V, const int E, const double bw_B, const double w_B,  const int d_num, const double dBWB  ){

  typedef float T;
  
  compressed_sparse_row_graph<int,int> graph;
  set<pair<int, int> > hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<T>  caps;
  vector<double> weights;

  std::default_random_engine generator;
  std::uniform_real_distribution<T> disBW(0.0,bw_B);
  
  std::uniform_real_distribution<T> disWS(0.0,w_B);

  std::uniform_real_distribution<T> disDBW(0.0, dBWB);


  int i = 0;
  int src, snk;
  // double weight,
  double bw;
  pair<int, int> temp;
  while (i < E) {
    src = rand() % V;
    snk = rand() % V;       
    while (src==snk) {
      snk = rand() % V;       
    }
    temp.first = src;
    temp.second = snk;
    if (hasSet.find(temp) == hasSet.end()) {
      i++;

      bw=disBW( generator );
      // weight = disWS(generator);

      hasSet.insert(temp);
      srcs.push_back(src);
      snks.push_back(snk);
      weights.push_back(1);
      caps.push_back( bw );

    }
  }
  i=0;

  typedef CG<compressed_sparse_row_graph<int, int>,  T> CG_T;
  vector<CG_T::Demand> demands;
  while (i< d_num) {
    src = rand() % V;
    snk = rand() % V;
    while(src==snk){
      snk = rand() % V;
    }
    bw=disDBW( generator );
    CG_T::Demand d;
    d.src=src;
    d.snk=snk;
    d.bandwidth=bw;
    demands.push_back( d );
    i++;
  }
  vector<int> ws( snks.size(  ), 1 );
  
  graph.initial( srcs, snks , ws);
  CG_T cg( graph,weights, caps, demands );
  cg.setInfo( 1 );
  cg.solve(  );
  
}

// void Case1(  ){
//   compressed_sparse_row_graph<float > graph;

//   vector<int> srcs;
//   vector<int> snks;
//   vector<float> weights;
//   srcs.push_back( 0 );
//   snks.push_back( 1 );

//   srcs.push_back( 0 );
//   snks.push_back( 2 );

//   srcs.push_back( 1 );
//   snks.push_back( 3 );

//   srcs.push_back( 2 );
//   snks.push_back( 3 );

//   srcs.push_back( 2 );
//   snks.push_back( 4 );

//   srcs.push_back( 1 );
//   snks.push_back( 4 );

//   srcs.push_back( 3 );
//   snks.push_back( 4 );

//   int num=7;
//   for (int i = 0; i < num; i++) {
//     weights.push_back( 1 );

//   }

//   graph.initial( srcs, snks, weights );
//   vector<int> path;
//   int re =graph.getShortPath( 0,4, path );
//   if( re>0 )
//     graph.printPath( path );

// }


void case2(void) {
  compressed_sparse_row_graph< float> graph(5);

  vector<int> srcs;
  vector<int> snks;
  vector<float> weights;

  srcs.push_back(0);
  snks.push_back(1);
  weights.push_back(3);

  srcs.push_back(1);
  snks.push_back(2);
  weights.push_back(4);

  srcs.push_back(2);
  snks.push_back(3);
  weights.push_back(1);

  srcs.push_back(4);
  snks.push_back(3);
  weights.push_back(2);

  srcs.push_back(5);
  snks.push_back(4);
  weights.push_back(3);

  srcs.push_back(0);
  snks.push_back(5);
  weights.push_back(2);

  srcs.push_back(5);
  snks.push_back(1);
  weights.push_back(1);

  srcs.push_back(5);
  snks.push_back(2);
  weights.push_back(2);

  srcs.push_back(2);
  snks.push_back(4);
  weights.push_back(2);

  graph.initial(srcs, snks, weights);
  graph.compute_allPair_shortest_path();
  vector<vector<int> > paths;
  graph.getShortPath(0, 3, paths);

  for (size_t i = 0; i < paths.size(); i++) {
    graph.printPath(paths[i]);
    std::cout << std::endl;
  }

  //  graph.increaseLinkWeight( 2, 3);
}

void randGraph(const int V, const int E, const double WW) {
  compressed_sparse_row_graph< double, int> graph;
  graph.setInfi(1000000000);

  set<pair<int, int> > hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<double> weights;

  int W= ( int )(WW+0.1);
  int i = 0;
  size_t src, snk;
  int weight;
  pair<int, int> temp;
  while (i < E) {
    src = rand() % V;
    snk = rand() % V;
    temp.first = src;
    temp.second = snk;
    if (hasSet.find(temp) == hasSet.end()) {
      i++;
      weight =  rand()%W+1 ;
      hasSet.insert(temp);
      srcs.push_back(src);
      snks.push_back(snk);
      weights.push_back(weight);
    }
  }

  graph.initial(srcs, snks, weights);
  // compressed_sparse_row_graph< double>::vertex_map<double> vmap=
  //     graph.get_vertex_map<double>(  );
  // vmap[ 1 ]=122.22;
  // vmap[ 2 ]=122.22;
  
  clock_t start=clock(  );
  // graph.compute_allPair_shortest_path();
  std::cout << ( (double)(clock(  )-start))/CLOCKS_PER_SEC<<"  seconds" << std::endl;
  size_t tlen=0;
  size_t tnum=0;
  vector<vector<int> > paths;
  

  for (i = 0; i < 1000; i++) {
    set<int> pset;
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }
    paths.clear();
    vector<int> path;
    if(graph.compute_shortest_path_dijkstra( src, snk, path )){
      if (!graph.isValidatePath(src, snk, path))
        std::cout << "error" << std::endl;
    }
    
    
    // graph.getShortPath(src, snk, paths);

    // tnum++;
    // for (size_t i = 0; i < paths.size(); i++) {
    //   pset.insert(paths[ i ].begin(  ), paths[ i ].end(  )  );

    //   if (!graph.isValidatePath(src, snk, paths[i]))
    //     std::cout << "error" << std::endl;

    // }
    // tlen+=pset.size(  );
  }
  std::cout <<tnum<<" pair "<< tlen<< " length"<<"  mean length: "<< (tlen/( tnum+0.01 )) << std::endl;

  // graph.increaseLinkWeight( V/2, 10 );
  // vector<int> path;
  // int re =graph.getShortPath( 0,V/2, path );
  // if( re>0 ){
  //   graph.printPath( path );
  // }
}

void randbiGraph( const int V, const int E, const double WW ){

  compressed_sparse_row_graph< double, int> graph;
  undir_graph< double, int> bgraph;
  graph.setInfi(100000);

  set<pair<int, int> > hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<double> weights;

  vector<int> bsrcs;
  vector<int> bsnks;
  vector<double> bweights;

  int W= ( int )(WW+0.1);
  int i = 0;
  int src, snk;
  int weight;
  pair<int, int> temp;
  while (i < E) {
    src = rand() % V;
    snk = rand() % V;       
    while (src==snk) {
      snk = rand() % V;       
    }
    temp.first = src;
    temp.second = snk;
    if (hasSet.find(temp) == hasSet.end()) {
      i++;
      weight =  rand()%W+1 ;
      hasSet.insert(temp);
      srcs.push_back(src);
      snks.push_back(snk);
      weights.push_back(weight);

      bsrcs.push_back(src);
      bsnks.push_back(snk);
      bweights.push_back(weight);

      bsrcs.push_back(snk);
      bsnks.push_back(src);
      bweights.push_back(weight);

    }
  }

  graph.initial(bsrcs, bsnks, bweights, false);
  bgraph.initial(srcs, snks  );
  
  clock_t start=clock(  );

  std::cout << ( (double)(clock(  )-start))/CLOCKS_PER_SEC<<"  seconds" << std::endl;
  size_t tlen=0;
  size_t tnum=0;
  vector<vector<int> > paths;
  int bsucc=0;
  int succ=0;
  double bsumdis=0;
  double sumdis=0;
  vector<pair<int, int> >task;
  vector<int> all_dis( 1000 );
  double start1=cpuTime(  )  ;
  for (i = 0; i < 1000; i++) {
    set<int> pset;
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }
    task.push_back( make_pair( src, snk ) );
 
    vector<int> path;
    if( bidijkstra_shortest_path( graph, bweights, src, snk, path , 10000 )){
      // if(graph.bicompute_shortest_path_dijkstra( src, snk, path )){
      all_dis[ i ]=graph.path_cost( path );
      bsumdis+=graph.path_cost( path );
      bsucc++;
      
      if (!graph.isValidatePath(src, snk, path))
        std::cout << "error" << std::endl;
    }

  }
  
  std::cout <<tnum<<" pair "<< tlen<< " length"<<"  mean length: "<< (tlen/( tnum+0.01 )) << std::endl;
  std::cout << "bi method success  "<<bsucc << std::endl;
  std::cout << "bi method time  "<<cpuTime(  )-start1 << std::endl;
  start1=cpuTime(  )  ;
  i=0;  
  for( vector<pair<int, int> > ::iterator it =task.begin( ); it!= task.end(  ); it++  ){
    vector<int> path, path1, path2, path3;
    src=it->first;
    snk=it->second;
    if( bidijkstra_shortest_path( bgraph, weights, src, snk, path , 10000 )){
      // graph.compute_shortest_path_dijkstra1( src, snk, path1 );
      // dijkstra_shortest_path( bgraph, weights, 10000, src, snk, path2  );
      // bidijkstra_shortest_path( graph, bweights, 10000, src, snk, path3  );

      // int temp=graph.path_cost( path );
      // // int temp1=graph.path_cost( path1 );
      
      int temp=path_cost(weights, path, 1 );
      // int temp3=graph.path_cost( path3);
      // assert( temp==temp1 );
      // assert( temp2==temp1 );
      // assert( temp2==temp3 );
      sumdis+=temp;
      succ++;
      // assert( graph.isValidatePath(src, snk, path) );
      // assert( graph.isValidatePath(src, snk, path1) );
      // if (!graph.isValidatePath(src, snk, path))
      //   std::cout << "error" << std::endl;
    }
    i++;
    
  }
  std::cout << " method success  "<<succ << std::endl;
  std::cout << bsumdis/sumdis << std::endl;
  std::cout << "bi method time  "<<cpuTime(  )-start1 << std::endl;
}

int main(int argc, char *argv[]) {
  
  // MCFexample2(  );
  randMCF( 60, 250, 300, 40, 20000, 3 );
  // //  Case1(  );
  // double start=cpuTime(  );
  // randbiGraph(20000, 100000, 20);
  
  // double end=cpuTime(  );
  // std::cout << "time "<<end-start<<" seconds" << "  "<<memUsedPeak( true )<< " M  "<<memUsed(  )<<"  M" << std::endl;
  // // case2(  );

  return 0;
}
