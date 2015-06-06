#include "composlist.h"
#include<set>
#include<vector>
#include <utility> 
#include<cstdlib>
#include<iostream>


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



void case2( void ){

  compressed_sparse_row_graph<size_t, float,float > graph;

  vector<size_t> srcs;
  vector<size_t> snks;
  vector<float> weights;

  srcs.push_back( 0 );
  snks.push_back( 1 );
  weights.push_back( 3 );



  srcs.push_back( 1);
  snks.push_back( 2 );
  weights.push_back( 4 );


  srcs.push_back( 2 );
  snks.push_back( 3 );
  weights.push_back( 1 );


  srcs.push_back( 4 );
  snks.push_back( 3 );
  weights.push_back( 2 );


  srcs.push_back( 5 );
  snks.push_back( 4 );
  weights.push_back( 3 );

  srcs.push_back( 0 );
  snks.push_back( 5 );
  weights.push_back( 2 );


  srcs.push_back( 5 );
  snks.push_back( 1 );
  weights.push_back( 1 );


  srcs.push_back( 5 );
  snks.push_back( 2 );
  weights.push_back( 2 );

  srcs.push_back( 2 );
  snks.push_back( 4 );
  weights.push_back( 2 );

  graph.initial( srcs, snks, weights );
  
  //  graph.increaseLinkWeight( 2, 3);
  
}

void randGraph( const int V , const  int E , const double W){


  compressed_sparse_row_graph<size_t, float,float > graph;

  set< pair<size_t, size_t> > hasSet;
  vector<size_t> srcs;
  vector<size_t> snks;
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
      weight=W*(rand( ))/RAND_MAX+0.1;
      hasSet.insert( temp );
      srcs.push_back( src );
      snks.push_back( snk );
      weights.push_back( weight );
      
    }

  }

  graph.initial( srcs, snks, weights );

  graph.compute_allPair_shortest_path(  );
  
  // for (i = 0; i < 1000; i++) {
  //   std::cout << i << std::endl;
  //   size_t link=rand( )% E;
  //   float f= rand( ) %100;
  //   graph.increaseLinkWeight( link, f);
  
  // }

  
  // graph.increaseLinkWeight( V/2, 10 );
  // vector<int> path;
  // int re =graph.getShortPath( 0,V/2, path );
  // if( re>0 ){
  //   graph.printPath( path );
  // }
  

}

int main(int argc, char *argv[])
{
  //  Case1(  );
  randGraph( 1000, 30000, 10 );

  
  return 0;
}
