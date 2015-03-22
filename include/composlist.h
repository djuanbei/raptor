/** 
 * only use in direct graph
 * 
 * 
 * @return 
 */
#ifndef COMPOSLIST_H
#define COMPOSLIST_H


#include<cstdio>
#include<vector>
#include<cassert>
#include<iostream>
using std::vector;


struct endElement{
  double weight;
  size_t snk;
  endElement(  ): weight( 0 ), snk( 0 ){
    
  }
};

struct startElement{
  size_t link;
  size_t src;
  startElement(  ):link( 0 ), src( 0 ){
    
  }
  
};

struct precedence
{
  size_t link;
  double weight;
  size_t vertex;
  precedence(  ):link( 0 ),  vertex( 0 ){
  }
  precedence  operator=( const precedence&  other){
    link=other.link ;
    weight=other.weight;
    vertex=other.vertex;
  }

};


class compressed_sparse_row_graph{

private:
  size_t vertex_num;
  size_t link_num;
  //out links
  vector<size_t> outIndex;
  vector< endElement> link_ends;
  //in links
  vector<size_t> inIndex;
  vector<startElement> link_starts;

  vector<size_t> link_map;// link_ends  --> orignal link
  vector<size_t> reLink_map; // orginal link --> link_end
  vector<precedence> shortPaths;  // allpair shortest path
  int _findSrc( const size_t link, size_t & src )const ;
  int _findSnk( const size_t link, size_t & snk ) const;

  inline int _getSrcSnk( const size_t link, size_t &src, size_t &snk  ) const{
    if(  0== _findSrc( link, src )) return 0;
    _findSrc( link, snk );
    return 1;
  }

  int _getShortPath( const size_t src, const size_t snk, vector<size_t> & path ) const;
  int _increaseLinkWeight(const  size_t link, const double inc );
  void compute_allPair_shortest_path(  );


public:
  compressed_sparse_row_graph( ):vertex_num( 0 ), link_num( 0 ){
  }

  void initial(vector<size_t> &srcs, vector< size_t > &snks, vector<double> &weights  );
  inline size_t getOutDegree(size_t vertex) const{
    assert( vertex< vertex_num );
    return outIndex[ vertex+1 ]-outIndex[ vertex ];
  }
  inline endElement & getLink(size_t vertex, size_t k ) {
    assert(  vertex< vertex_num && k< getOutDegree( vertex )  );
    return  link_ends[outIndex[vertex ]+k ];
  }
  inline int findSrc( const size_t link, size_t & src )const {
      return _findSrc( reLink_map[ link ], src );
  }
  inline int findSnk( const size_t link, size_t & snk ) const{
    return _findSnk(reLink_map[ link ], snk  );
  }

  inline int getSrcSnk( const size_t link, size_t &src, size_t &snk  ) const{
    return _getSrcSnk(reLink_map[link], src, snk  );
    
  }

  inline int increaseLinkWeight(const  size_t link, const double inc ){
    return _increaseLinkWeight(reLink_map[link], inc ) ;

  }

  int  getShortPath( const size_t src, const size_t snk, vector<size_t> & path ) const ;

  static void printPath( vector<size_t>&path ){
    if(0== path.size(  ) ) return;
    size_t i=0;

    std::cout<<path[ 0 ]<<std::endl;
    i++;
    while( i<path.size(  ) ){
      std::cout<<"by " << path[ i ]<<std::endl;
      i++;
      std::cout<<path[ i ]<<std::endl;
      i++;
    }
  }
  
};


#endif
