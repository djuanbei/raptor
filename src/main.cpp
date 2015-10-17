#include "graph.h"
#include "System.h"
#include <set>
#include <vector>
#include <utility>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace fast_graph;

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

int main(int argc, char *argv[]) {
  //  Case1(  );
  double start=cpuTime(  );
  randGraph(100000, 300000, 10);
  
  double end=cpuTime(  );
  std::cout << "time "<<end-start<<" seconds" << "  "<<memUsedPeak( true )<< " M  "<<memUsed(  )<<"  M" << std::endl;
  // case2(  );

  return 0;
}
