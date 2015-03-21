#ifndef FAST_GRAPH_H
#defind FAST_GRAPH_H
#include<cstdlib>

struct undirectedS {
  static bool is_directed() { return false; }
  static bool is_bidirectional() { return false; }
  static const unsigned long direction_type = 0;
};

struct directedS {
  static bool is_directed() { return true; }
  static bool is_bidirectional() { return false; }
  static const unsigned long direction_type = 1;
};



template <typename DIRECTION = directedS    >
class compressed_sparse_row_graph {
public:
  typedef unsigned long size_type;  
  typedef double value_type;
  typedef double cap_type;
  compressed_sparse_row_graph( ):m( 0 ),n( 0 ),index( NULL ),edge_points( NULL ),pair_edge_loc( NULL ),edge_orignal_id_map( NULL ){
  }
  ~compressed_sparse_row_graph( ) {
    clear( ) ;
  }
  void clear( ) {
    if(NULL!= index ){
      free( index );
      index=NULL;
    }
    if( NULL!=edge_points ){
      free( edge_points );
      edge_points=NULL;
  
    }
    if( NULL!=pair_edge_loc ){
      free( pair_edge_loc );
      pair_edge_loc=NULL;
    }
    
    if(NULL!=edge_orignal_id_map ){
      free( edge_orignal_id_map );
      edge_orignal_id_map=NULL;
    }
    m=n=0;
  }


  void init( size_type vertex_num, size_type edge_num, size_type * srcs, size_type *dests, value_type * weights, cap_type *caps  ){
    n=vertex_num;
    m=edge_num;
    size_type a_size=DIRECTED::is_directed( )? m: 2*m ;
  }

private:
  struct element
  {
    size_type sink;
    value_type weight;
    cap_type cap;
  };
  size_type m; // of edges
  size_type n;// of vertices
  compressed_sparse_row_graph & operator =(const compressed_sparse_row_graph& rhs  ){
    
  }

  // Indexes into the edge-arrays for the adjacencies.  Its size is |V|+1.  A
  // vertex, v, has neighbors in {src,end}_points from {src,end}_points[v] to
  // {src,end}_points[v+1].
  size_type* index;

  // edge_points is edges rearranged  from the vertices order
  element *edge_points;
  size_type *pair_edge_loc;

  size_type *edge_orignal_id_map;
    
}




#endif
