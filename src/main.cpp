#include "composlist.h"
#include<set>
#include<vector>
#include <utility> 
#include<cstdlib>
#include<iostream>

using std::set;
using std::vector;
using std::pair;



void Case1(  ){
  compressed_sparse_row_graph graph;

  vector<size_t> srcs;
  vector<size_t> snks;
  vector<double> weights;
  srcs.push_back( 0 );
  snks.push_back( 1 );

  srcs.push_back( 0 );
  snks.push_back( 2 );
    
  srcs.push_back( 1 );
  snks.push_back( 3 );
    
  srcs.push_back( 2 );
  snks.push_back( 3 );
    
  srcs.push_back( 2 );
  snks.push_back( 4 );
    
  srcs.push_back( 1 );
  snks.push_back( 4 );
    
  srcs.push_back( 3 );
  snks.push_back( 4 );
    
  int num=7;
  for (int i = 0; i < num; i++) {
    weights.push_back( 1 );
  
  }

  graph.initial( srcs, snks, weights );
  vector<size_t> path;
  int re =graph.getShortPath( 0,4, path );
  graph.printPath( path );

}



void randGraph( const int N , const  int V , const double W){
  compressed_sparse_row_graph graph;

  set< pair<size_t, size_t> > hasSet;
  vector<size_t> srcs;
  vector<size_t> snks;
  vector<double> weights;

  int i=0;
  size_t src, snk;
  double weight;
  pair<size_t, size_t> temp;
  while( i< V ){
    src=rand( )% N;
    snk=rand(  ) %N;
    temp.first=src;
    temp.second=snk;
    if( hasSet.find( temp )==hasSet.end(  ) ){
      i++;
      weight=W*(rand( )+0.1)/RAND_MAX;
      hasSet.insert( temp );
      srcs.push_back( src );
      snks.push_back( snk );
      weights.push_back( weight );
      
    }

  }
  graph.initial( srcs, snks, weights );
  vector<size_t> path;
  int re =graph.getShortPath( 0,N/2, path );
  if( re>0 ){
    graph.printPath( path );
  }
}

int main(int argc, char *argv[])
{
  //  Case1(  );
  randGraph( 500, 100000, 100 );
  
  return 0;
}
