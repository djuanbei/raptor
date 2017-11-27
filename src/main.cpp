#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <utility>
#include <vector>
#include "System.h"
#include "config.h"
#include "graph.h"
#include "graph.hpp"
#include "graphalg.hpp"

#include "mcfcg.hpp"

#include "csvReader.hpp"
#include "sparse.h"

using namespace std;
using namespace raptor;
using namespace mcmcf;
using namespace sparse;

void testSparse() {
  vector<sparseMatrixElem> elements;
  sparseMatrixElem temp;
  temp.row = 0;
  temp.column = 1;
  temp.value = 1;
  elements.push_back(temp);

  temp.row = 1;
  temp.column = 0;
  elements.push_back(temp);

  temp.row = 2;
  temp.column = 2;
  elements.push_back(temp);
  SparseSolver solver;

  solver.update(elements);
  double b[3];
  b[0] = 1;
  b[1] = 0;
  b[2] = 0;

  solver.locSolver(b);
}

void example2(void) {
  vector<int> srcs;
  vector<int> snks;
  vector<double> weights, caps;
  srcs.push_back(0);
  snks.push_back(1);
  weights.push_back(2);
  caps.push_back(1);

  srcs.push_back(0);
  snks.push_back(1);
  weights.push_back(1);
  caps.push_back(2);

  srcs.push_back(2);
  snks.push_back(0);
  weights.push_back(1);
  caps.push_back(2);

  srcs.push_back(2);
  snks.push_back(1);
  weights.push_back(4);
  caps.push_back(2);

  simple_graph graph;

  typedef CG<simple_graph, double> CG_T;

  vector<Demand<double>> demands;
  Demand<double> d1;
  d1.src = 0;
  d1.snk = 1;
  d1.bandwidth = 1;
  demands.push_back(d1);
  Demand<double> d2;
  d2.src = 2;
  d2.snk = 1;
  d2.bandwidth = 2;

  demands.push_back(d2);

  graph.initial(srcs, snks);
  CG_T cg(graph, weights, caps, demands);
  cg.setInfo(1);
  cg.solve();
}

void randMCF(int solver, const int V, const int E, const double bw_B,
             const double w_B, const int d_num, const double dBWB) {
  typedef double T;

  simple_graph graph;
  set<pair<int, int>> hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<T> caps;
  vector<double> weights;

  std::default_random_engine generator;
  std::uniform_real_distribution<T> disBW(0.0, bw_B);

  std::uniform_real_distribution<T> disWS(0.0, w_B);

  std::uniform_real_distribution<T> disDBW(0.0, dBWB);

  int i = 0;
  int src, snk;
  double weight;
  double bw;
  pair<int, int> temp;
  while (i < E) {
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }
    temp.first = src;
    temp.second = snk;
    if (hasSet.find(temp) == hasSet.end()) {
      i++;

      bw = (int)disBW(generator) + 1;
      weight = (int)disWS(generator) + 1;

      hasSet.insert(temp);
      srcs.push_back(src);
      snks.push_back(snk);
      weights.push_back(weight);
      caps.push_back(bw);
    }
  }
  i = 0;

  typedef CG<simple_graph, T> CG_T;
  vector<Demand<T>> demands;
  while (i < d_num) {
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }
    bw = (int)disDBW(generator) + 1;
    Demand<T> d;
    d.src = src;
    d.snk = snk;
    d.bandwidth = bw;
    demands.push_back(d);
    i++;
  }
  vector<int> ws(snks.size(), 1);

  graph.initial(srcs, snks);
  CG_T cg(graph, weights, caps, demands);
  cg.setInfo(2);
  if (solver == 0) {
    cg.setLUSOLVER(KLU);
  } else {
    cg.setLUSOLVER(LAPACK);
  }
  // cg.writeKsptoCNF(10, "test.cnf");
  cg.solve();
}

// void Case1(  ){
//   directed_graph<float > graph;

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

void randGraph(const int V, const int E, const double WW) {
  directed_graph<double, int> graph;
  graph.setInfi(1000000000);

  set<pair<int, int>> hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<double> weights;

  int W = (int)(WW + 0.1);
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
      weight = rand() % W + 1;
      hasSet.insert(temp);
      srcs.push_back(src);
      snks.push_back(snk);
      weights.push_back(weight);
    }
  }

  graph.initial(srcs, snks, weights);
  // directed_graph< double>::vertex_map<double> vmap=
  //     graph.get_vertex_map<double>(  );
  // vmap[ 1 ]=122.22;
  // vmap[ 2 ]=122.22;

  size_t tlen = 0;
  size_t tnum = 0;
  double slen = 0;
  size_t num = 10000;

  vector<vector<int>> paths;
  vector<int> dsrcs;
  snk = rand() % V;
  vector<double> dis(graph.getVertex_num());
  vector<int> preLink;

  dijkstra_shortest_retree(graph, weights, snk, preLink, dis, 100000000.0);

  clock_t start = clock();
  // graph.compute_allPair_shortest_path();
  std::cout << ((double)(clock() - start)) / CLOCKS_PER_SEC << "  seconds"
            << std::endl;

  for (size_t i = 0; i < num; i++) {
    set<int> pset;
    src = rand() % V;

    while (src == snk) {
      src = rand() % V;
    }
    dsrcs.push_back(src);
    paths.clear();
    vector<int> path;
    if (bidijkstra_shortest_path(graph, weights, src, snk, path, 1000000.0)) {
      if (!graph.isValidatePath(src, snk, path))
        std::cout << "error" << std::endl;
      tlen += path.size();
      tnum++;
      slen += path_cost(weights, path, 0.0);
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
  double d = clock() - start;
  std::cout << "time: " << d << std::endl;

  std::cout << tnum << " pair " << tlen << " length"
            << "  mean length: " << (tlen / (tnum + 0.01))
            << "        :    sum length:  " << slen << std::endl;
  tlen = 0;
  tnum = 0;
  slen = 0;
  start = clock();
  for (size_t i = 0; i < num; i++) {
    src = dsrcs[i];
    paths.clear();
    vector<int> path;

    if (astar_shortest_path(graph, weights, dis, src, snk, path, 100000.0)) {
      if (!graph.isValidatePath(src, snk, path))
        std::cout << "error" << std::endl;
      tlen += path.size();
      tnum++;
      slen += path_cost(weights, path, 0.0);
    }
  }

  d = clock() - start;
  std::cout << "time: " << d << std::endl;

  std::cout << tnum << " pair " << tlen << " length"
            << "  mean length: " << (tlen / (tnum + 0.01))
            << "        :    sum length:  " << slen << std::endl;
  // graph.increaseLinkWeight( V/2, 10 );
  // vector<int> path;
  // int re =graph.getShortPath( 0,V/2, path );
  // if( re>0 ){
  //   graph.printPath( path );
  // }
}

void randbiGraph(const int V, const int E, const double WW) {
  directed_graph<double, int> graph;
  undirected_graph<double, int> bgraph;
  graph.setInfi(100000);

  set<pair<int, int>> hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<double> weights;

  vector<int> bsrcs;
  vector<int> bsnks;
  vector<double> bweights;

  int W = (int)(WW + 0.1);
  int i = 0;
  int src, snk;
  int weight;
  pair<int, int> temp;
  while (i < E) {
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }
    temp.first = src;
    temp.second = snk;
    if (hasSet.find(temp) == hasSet.end()) {
      i++;
      weight = rand() % W + 1;
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
  bgraph.initial(srcs, snks);

  clock_t start = clock();

  std::cout << ((double)(clock() - start)) / CLOCKS_PER_SEC << "  seconds"
            << std::endl;
  size_t tlen = 0;
  size_t tnum = 0;
  vector<vector<int>> paths;
  int bsucc = 0;
  int succ = 0;
  double bsumdis = 0;
  double sumdis = 0;
  vector<pair<int, int>> task;
  vector<int> all_dis(1000);
  double start1 = cpuTime();
  for (i = 0; i < 1000; i++) {
    set<int> pset;
    src = rand() % V;
    snk = rand() % V;
    while (src == snk) {
      snk = rand() % V;
    }
    task.push_back(make_pair(src, snk));

    vector<int> path;
    if (bidijkstra_shortest_path(graph, bweights, src, snk, path, 10000)) {
      // if(graph.bicompute_shortest_path_dijkstra( src, snk, path )){
      all_dis[i] = graph.path_cost(path);
      bsumdis += graph.path_cost(path);
      bsucc++;

      if (!graph.isValidatePath(src, snk, path))
        std::cout << "error" << std::endl;
    }
  }

  std::cout << tnum << " pair " << tlen << " length"
            << "  mean length: " << (tlen / (tnum + 0.01)) << std::endl;
  std::cout << "bi method success  " << bsucc << std::endl;
  std::cout << "bi method time  " << cpuTime() - start1 << std::endl;
  start1 = cpuTime();
  i = 0;
  for (vector<pair<int, int>>::iterator it = task.begin(); it != task.end();
       it++) {
    vector<int> path, path1, path2, path3;
    src = it->first;
    snk = it->second;
    if (bidijkstra_shortest_path(bgraph, weights, src, snk, path, 10000)) {
      // graph.compute_shortest_path_dijkstra1( src, snk, path1 );
      // dijkstra_shortest_path( bgraph, weights, 10000, src, snk, path2
      // );
      // bidijkstra_shortest_path( graph, bweights, 10000, src, snk, path3
      // );

      // int temp=graph.path_cost( path );
      // // int temp1=graph.path_cost( path1 );

      int temp = path_cost(weights, path, 1);
      // int temp3=graph.path_cost( path3);
      // assert( temp==temp1 );
      // assert( temp2==temp1 );
      // assert( temp2==temp3 );
      sumdis += temp;
      succ++;
      // assert( graph.isValidatePath(src, snk, path) );
      // assert( graph.isValidatePath(src, snk, path1) );
      // if (!graph.isValidatePath(src, snk, path))
      //   std::cout << "error" << std::endl;
    }
    i++;
  }
  std::cout << " method success  " << succ << std::endl;
  std::cout << bsumdis / sumdis << std::endl;
  std::cout << "bi method time  " << cpuTime() - start1 << std::endl;
}

void testAns(char *filename) {
  ifstream ifs;
  ifs.open(filename, std::ifstream::in);
  csv_istream csv_in(ifs);
  vector<int> srcs, snks;
  vector<double> weights;
  int src, snk;
  double w;
  while (csv_in >> src && csv_in >> snk) {
    csv_in >> w;
    srcs.push_back(src);
    snks.push_back(snk);
    srcs.push_back(snk);
    snks.push_back(src);
    weights.push_back(w);
    weights.push_back(w);
  }
  double infi_value = 1000000000;
  directed_graph<double, int> graph;
  graph.initial(srcs, snks, weights);
  inc_ksp::yen_ksp<directed_graph<double, int>, vector<double>, double> yen(
      graph, weights, infi_value);

  inc_ksp::yen_next_path<directed_graph<double, int>, vector<double>, double>
      next_p = yen.next_path(15, 0);

  int path_num = 0;
  vector<int> path;
  set<vector<int>> paths;
  while (next_p.next_path(path)) {
    if (!isSimplePath(graph, 15, 0, path)) {
      cout << "error" << endl;
    }
    sort(path.begin(), path.end());
    if (paths.find(path) != paths.end()) {
      cout << "error" << endl;
    }
    paths.insert(path);
    path_num++;
  }
  cout << "the number of simple path connect Hawii and Hartford is: "
       << path_num << endl;
}

int main(int argc, char *argv[]) {
  testSparse();
  return 0;
  // randMCF(0, 20, 100, 400, 10, 60, 100);
  // example2();
  // //  testAns(argv[1]);
  // return 0;

  // ifstream ifs;
  // ifs.open(argv[1], std::ifstream::in);
  // csv_istream csv_in(ifs);

  // string temp;
  // while (csv_in >> temp) {
  //   cout << temp << endl;
  // }

  // ifs.close();

  // return 0;

  // randGraph(20000, 100000, 20);
  // return 0;

  // MCFexample2(  );
  // randMCF(4, 8, 200, 40, 6, 50);
  cout << "case, using time(s), success rat, object value, iteration number, "
          "empty iteration,  computing shortest path use time(s), solving "
          "linear equation solve use time(s)"
       << endl;
  if (argc > 1) {
    for (int i = 1; i < 40; i += 2) {
      // cout<<"************************************"<<endl;
      // cout<<1000*i<<endl;
      // cout<<"************************************"<<endl;
      cout << "case_" << 1000 * i << "_" << 5000 * i << ",";
      randMCF(1, 1000 * i, 5000 * i, 300, 10, 1000, 100);
    }
  } else {
    for (int i = 1; i < 40; i += 2) {
      // cout<<"************************************"<<endl;
      // cout<<1000*i<<endl;
      // cout<<"************************************"<<endl;
      cout << "case_" << 1000 * i << "_" << 5000 * i << ",";
      randMCF(0, 1000 * i, 5000 * i, 300, 10, 1000, 100);
    }
  }
  // //  Case1(  );
  // double start=cpuTime(  );
  // randbiGraph(20000, 100000, 20);

  // double end=cpuTime(  );
  // std::cout << "time "<<end-start<<" seconds" << "  "<<memUsedPeak( true
  // )<< " M  "<<memUsed(  )<<"  M" << std::endl;
  // // case2(  );

  return 0;
}
