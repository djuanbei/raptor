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
  tempElment dummy;

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
    dummy.id=i;
    dummy.src=srcs[ i ];
    tContian.push_back( dummy );
  }
  std::sort( tContian.begin(  ), tContian.end(  ), tempElmentCmp  );

  index.push_back( 0 );
  i=0;
  while(i< tContian[ 0 ].src  ){
    index.push_back( 0 );
    i++;
  }

  endElement temp;
  temp.snk=snks[tContian[ 0 ].id ];
  temp.weight=weights[tContian[ 0 ].id ];

  link_ends.push_back( temp );
  link_map[ 0 ]=tContian[ 0 ].id;
  reLink_map[ tContian[ 0 ].id ]=0;
  
  for( i=1; i< tContian.size( ) ;i++ ){
    if(  tContian[ i ].src!= tContian[ i-1 ].src ){
      for(j=tContian[ i-1 ].src; j< tContian[ i ].src; j++ ){
        index.push_back(link_map.size( ) );
      }
    }
    temp.snk=snks[ tContian[ i ].id ];
    temp.weight=weights[ tContian[ i ].id ];

    link_ends.push_back( temp );
    link_map[ i ]=tContian[ i ].id;
    reLink_map[tContian[ i ].id  ]=i;
  }

  index.push_back( link_map.size( ) );

  compute_allPair_shortest_path(  );
}

void compressed_sparse_row_graph::compute_allPair_shortest_path(  ){

  size_t i;
  precedence infPre;
  infPre.link=link_num+1;
  infPre.weight=numeric_limits<double>::max(  );
  infPre.vertex=vertex_num+1;

  shortPaths.resize(vertex_num*vertex_num, infPre );
#pragma omp parallel for
  for( i=0; i< vertex_num; i++ ){
    const  size_t shift=i*vertex_num;
    size_t j, current, out;
    double w;

    shortPaths[shift+i].weight=0;
    queue<size_t>  wait;
    wait.push(i );
    while( !wait.empty(  ) ){
      current=wait.front(  );
      wait.pop(  );
      out=getOutDegree( current );
      for( j=0; j< out; j++ ){

        precedence &dummy=shortPaths[shift + link_ends[index[current]+j].snk];
        w=shortPaths[ shift+current ].weight+ link_ends[index[ current ]+j].weight;

        if( w <dummy.weight   ){
          dummy.link= index[current]+j;
          dummy.weight =w;
          dummy.vertex=current;
          wait.push( link_ends[index[ current]+j].snk );

        }
      }
    }
  }

}

int compressed_sparse_row_graph::increaseLinkWeight(const  size_t i, const double inc ){
  
  return 0;
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
vector<size_t> & getShortPath( const size_t src, const size_t snk ){
  assert( src< vertex_num && snk < vertex_num );
  assert( src!=snk );
  vector<size_t> path;
  const size_t shift= src*vertex_num;
  if(shortPaths[ shift +snk ].vertex> vertex_num )
    return  path;

  
  size_t  current =snk;
  while( current!=src ){
    path.push_back( current );
    path.push_back(shortPaths[ shift +current ].link  );
    current=shortPaths[ shift +current ].vertex;
  }
  path.push_back( current );
  std::reverse(path.begin(  ), path.end(  )  );

  return path;
  
}
