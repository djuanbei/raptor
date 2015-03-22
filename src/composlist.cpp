#include<algorithm>
#include<limits> 
#include<queue>
#include"composlist.h"

using std::numeric_limits;
using std::queue;

struct tempElment
{
  size_t id;
  size_t src;
};

static bool tempElmentCmp(const tempElment &lhs, const tempElment &rhs  ) {
  return lhs.src<rhs.src;
}



void compressed_sparse_row_graph::initial(vector<size_t> &srcs, vector< size_t > &snks, vector<double> &weights   ){

  assert(srcs.size( )==  snks.size(  ) && srcs.size(  )== weights.size(  ) );
  if ( 0== srcs.size(  ) ) return ;
  link_num=srcs.size(  );
  link_map.resize(link_num  );
  reLink_map.resize( link_num );

  vector<tempElment> tContian;
  size_t i ,j, ver_max;
  tempElment dummy_pred;

  ver_max=srcs[ 0 ];

  for( i=0; i< link_num; i++ ){
    if( ver_max< srcs[ i ] ){
      ver_max=srcs[ i ];      
    }
    if( ver_max< snks[ i ] ){
      ver_max=snks[ i ];
    }
  }
  vertex_num=ver_max+1;

  for(i=0 ; i< srcs.size( ) ; i++ ){
    dummy_pred.id=i;
    dummy_pred.src=srcs[ i ];
    tContian.push_back( dummy_pred );
  }
  std::sort( tContian.begin(  ), tContian.end(  ), tempElmentCmp  );

  index.push_back( 0 );
  i=0;
  while(i< tContian[ 0 ].src  ){
    index.push_back( 0 );
    i++;
  }

  endElement temp;
  temp.weight=weights[tContian[ 0 ].id ];
  temp.snk=snks[tContian[ 0 ].id ];


  link_ends.push_back( temp );
  link_map[ 0 ]=tContian[ 0 ].id;
  reLink_map[ tContian[ 0 ].id ]=0;
  
  for( i=1; i< tContian.size( ) ;i++ ){
    if(  tContian[ i ].src!= tContian[ i-1 ].src ){
      for(j=tContian[ i-1 ].src; j< tContian[ i ].src; j++ ){
        index.push_back(link_ends.size( ) );
      }
    }
    temp.weight=weights[ tContian[ i ].id ];
    temp.snk=snks[ tContian[ i ].id ];


    link_ends.push_back( temp );
    link_map[ i ]=tContian[ i ].id;
    reLink_map[tContian[ i ].id  ]=i;
  }

  for( j=tContian[ i-1 ].src; j< vertex_num; j++ ){
      index.push_back(link_ends.size( ) );
  }


  compute_allPair_shortest_path(  );
}

void compressed_sparse_row_graph::compute_allPair_shortest_path(  ){

  size_t i;
  precedence infPre;
  infPre.link=link_num+1;
  infPre.weight=numeric_limits<double>::max( )/10;
  infPre.vertex=vertex_num+10;

  shortPaths.resize(vertex_num*vertex_num, infPre );
#pragma omp parallel for
  for( i=0; i< vertex_num; i++ ){
    const  size_t shift=i*vertex_num;
    size_t j, current, outDegree ;
  
    double weight;

    shortPaths[shift+i].weight=0;
    queue<size_t>  wait;
    wait.push(i );
    while( !wait.empty(  ) ){
      current=wait.front(  );
      wait.pop(  );
      outDegree=getOutDegree( current );
      for( j=0; j< outDegree; j++ ){
        endElement &neighbour=link_ends[index[current]+j];
        precedence &dummy_pred=shortPaths[shift + neighbour.snk];
        weight=shortPaths[ shift+current ].weight+ neighbour.weight;

        if( weight <dummy_pred.weight   ){
          dummy_pred.link= index[current]+j;
          dummy_pred.weight =weight;
          dummy_pred.vertex=current;
          wait.push( neighbour.snk );
        }
      }
    }
  }

}

int compressed_sparse_row_graph::increaseLinkWeight(const  size_t i, const double inc ){
  
  return 0;
}
int compressed_sparse_row_graph:: _getShortPath( const size_t src, const size_t snk, vector<size_t> & path ) const{
  assert( src< vertex_num && snk < vertex_num );
  assert( src!=snk );

  const size_t shift= src*vertex_num;
  if(shortPaths[ shift +snk ].vertex> vertex_num )
    return  0;
  
  size_t  current =snk;
  while( current!=src ){
    path.push_back( current );
    path.push_back( shortPaths[ shift +current ].link   );
    current=shortPaths[ shift +current ].vertex;
  }
  path.push_back( current );
  std::reverse(path.begin(  ), path.end(  )  );
  return 1;


}

/** 
 * 
 * 
 * @param src 
 * @param snk 
 * 
 * @return  the shortest path start from src and end with snk 
 * if there is no path connect with src and snk then  return an 
 * empty vector
 */
int compressed_sparse_row_graph::getShortPath( const size_t src, const size_t snk, vector<size_t> &path ) const{
  assert( src< vertex_num && snk < vertex_num );
  assert( src!=snk );
 int re=  _getShortPath( src, snk, path);
 if( re>0 ){
   size_t i=0;
   while( i+1< path.size(  ) ){
     path[ i+1 ]=link_map[ path[ i+1 ] ];
     i+=2;
   }
 }

 return re;
  
}
