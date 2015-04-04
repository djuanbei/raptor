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
  compressed_sparse_row_graph<int,float, float> graph;

  vector<int> srcs;
  vector<int> snks;
  vector<float> weights;
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
  vector<int> path;
  int re =graph.getShortPath( 0,4, path );
  if( re>0 )
    graph.printPath( path );

}


void randGraph( const int V , const  int E , const double W){
  compressed_sparse_row_graph<int, float,float > graph;
  graph.setInfi( 10000000 );

  set< pair<size_t, size_t> > hasSet;
  vector<int> srcs;
  vector<int> snks;
  vector<float> weights;

  int i=0;
  size_t src, snk;
  double weight;
  pair<size_t, size_t> temp;
  while( i< E ){
    src=rand( )% V;
    snk=rand( ) %V;
    temp.first=src;
    temp.second=snk;
    if( hasSet.find( temp )==hasSet.end(  ) ){
      i++;
      weight=W*(rand( ))/RAND_MAX;
      hasSet.insert( temp );
      srcs.push_back( src );
      snks.push_back( snk );
      weights.push_back( weight );
      
    }

  }

  graph.initial( srcs, snks, weights );
  // for (i = 0; i < 100; i++) {
  //   size_t link=rand( )% E;
  //   graph.increaseLinkWeight( link, 10 );
  
  // }

  
  // graph.increaseLinkWeight( V/2, 10 );
  // vector<size_t> path;
  // int re =graph.getShortPath( 0,V/2, path );
  // if( re>0 ){
  //   graph.printPath( path );
  // }

}

int main(int argc, char *argv[])
{
  //  Case1(  );
  randGraph( 1000, 100000, 10 );
  
  return 0;
}
