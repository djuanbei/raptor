/**
 * @file   heap.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Mon Apr  6 20:54:06 2015
 * 
 * @brief  the min heap used in dijstra shortest path algorithm with the graph only with nonnegtative 
 * edge weight
 * 
 * 
 */

#ifndef _HEAP_H
#define _HEAP_H
#include<vector>
#include <utility>      // std::pair, std::make_pair
using std::vector;
using std::pair;


template<class W, class E, class Comp >

class Heap{
private:
  vector<pair<W, E> > data; // Heap of Keys
  vector<int> indices;  // Each Key's position (index) in the Heap
  Comp                  lt;       // The data is a minimum-data with respect to this comparator
  int size;

  // Index "traversal" functions
  static inline int left  (int i) { return i*2+1; }
  static inline int right (int i) { return (i+1)<<1; }
  static inline int parent(int i) { return (i-1) >> 1; }
  void percolateUp(int i)
  {
    pair<W, E>   x  = data[i];
    int p  = parent(i);
        
    while (i != 0 && lt(x, data[p])){
      data[i]          = data[p];
      indices[data[i].second] = i;
      i                = p;
      p                = parent(p);
    }
    data   [i] = x;
    indices[x.second] = i;

  }

  void percolateDown(int i)
  {
    pair<W,E> x = data[i];
    while (left(i) <size){
      int child = right(i) < size && lt(data[right(i)], data[left(i)]) ? right(i) : left(i);
      if (!lt(data[child], x)) break;
      data[i]          = data[child];
      indices[data[i].second] = i;
      i                = child;
    }
    data   [i] = x;
    indices[x.second] = i;

  }
  
public:
  Heap(const Comp& c ) :  lt(c), size( 0 )  {}
  void setCap( int cap ){
    data.resize(cap );
    indices.resize( cap,-1 );
  }

  bool empty     ()          const { return 0==size; }
  void push(pair<W,E> k)
  {
    if(-1==indices[ k.second ]  ){
      data[ size ]=k;
      percolateUp(size);
      indices[ k.second ]=size;
      size++;
    }else{
      percolateUp(indices[ k.second ]);
    }
    

  }
  pair<W,E> & top(  ){
    return data[ 0 ];
  }

  void pop(  ){
    indices[data[0].second  ]=-1;
    data[0]          = data[ size-1 ];
    size--;
    if (size > 1) percolateDown(0);
 
  }

};


#endif
