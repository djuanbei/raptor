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
using std::vector;


struct endElement{
  double weight;
  size_t snk;
  endElement(  ): weight( 0 ), snk( 0 ){
    
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
  vector<size_t> index;
  vector< endElement> link_ends;
  vector<size_t> link_map;// link_ends  --> orignal link
  vector<size_t> reLink_map; // orginal link --> link_end
  vector<precedence> shortPaths;  // allpair shortest path
  void compute_allPair_shortest_path(  );


 public:
  compressed_sparse_row_graph( ):vertex_num( 0 ), link_num( 0 ){
  }

  void initial(vector<size_t> &srcs, vector< size_t > &snks, vector<double> &weights  );
  inline size_t getOutDegree(size_t i) const{
    assert( i< vertex_num );

    return index[ i+1 ]-index[ i ];
  }
  inline endElement & getLink(size_t i, size_t k ) {
    assert(  i< vertex_num && k< getOutDegree( i )  );
    return  link_ends[index[ i ]+k ];
  }
  int increaseLinkWeight(const  size_t i, const double inc );

  vector<size_t> & getShortPath( const size_t src, const size_t snk );
  
};


#endif
