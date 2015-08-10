#include "graph.h"
#include <set>
#include <vector>
#include <utility>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace fast_graph;

// void Case1(  ){
//   compressed_sparse_row_graph<int,float, float> graph;

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
  compressed_sparse_row_graph<size_t, float, float> graph(5);

  vector<size_t> srcs;
  vector<size_t> snks;
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
  vector<vector<size_t> > paths;
  graph.getShortPath(0, 3, paths);

  for (size_t i = 0; i < paths.size(); i++) {
    graph.printPath(paths[i]);
    std::cout << std::endl;
  }

  //  graph.increaseLinkWeight( 2, 3);
}

void randGraph(const int V, const int E, const double WW) {
  compressed_sparse_row_graph<size_t, int, int> graph(10);
  graph.setInfi(1000000000);

  set<pair<size_t, size_t> > hasSet;
  vector<size_t> srcs;
  vector<size_t> snks;
  vector<int> weights;

  int W= ( int )(WW+0.1);
  int i = 0;
  size_t src, snk;
  int weight;
  pair<size_t, size_t> temp;
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

  graph.compute_allPair_shortest_path();
  size_t tlen=0;
  size_t tnum=0;
  vector<vector<size_t> > paths;
  for (i = 0; i < 1000; i++) {
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }

    paths.clear();
    graph.getShortPath(src, snk, paths);

    tnum+=paths.size(  );
    for (size_t i = 0; i < paths.size(); i++) {
      tlen+=paths[i].size(  );
      if (!graph.isValidatePath(src, snk, paths[i]))
        std::cout << "error" << std::endl;
      // graph.printPath( paths[ i ] );
      //       std::cout << std::endl;
    }
  }
  std::cout <<tnum<<" pair "<< tlen<< " length"<<"  mean length: "<< (tlen/( tnum+0.0 )) << std::endl;

  // graph.increaseLinkWeight( V/2, 10 );
  // vector<int> path;
  // int re =graph.getShortPath( 0,V/2, path );
  // if( re>0 ){
  //   graph.printPath( path );
  // }
}

int main(int argc, char *argv[]) {
  //  Case1(  );
  randGraph(2000, 5000, 5);
  // case2(  );

  return 0;
}
