#include "composlist.h"
#include<set>
#include<vector>
#include <utility> 
#include<cstdlib>


using std::set;
using std::vector;
using std::pair;



void randGraph( const int N , const  int V , const double W){
  compressed_sparse_row_graph graph;

  set< pair<size_t, size_t> > hasSet;
  vector<size_t> srcs;
  vector<size_t> snks;
  vector<double> weights;

  int i=0;
  size_t src, snk;
  double w;
  pair<size_t, size_t> temp;
  while( i< V ){
    src=rand( )% N;
    snk=rand(  ) %N;
    temp.first=src;
    temp.second=snk;
    if( hasSet.find( temp )==hasSet.end(  ) ){
      i++;
      w=W*(rand( )+0.1)/RAND_MAX;
      hasSet.insert( temp );
      srcs.push_back( src );
      snks.push_back( snk );
      weights.push_back( w );
      
    }

  }
  graph.initial( srcs, snks, weights );

}

int main(int argc, char *argv[])
{
  randGraph( 100,1000, 10 );
  
  return 0;
}
