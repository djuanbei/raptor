
/**
 * @file   composlist.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Thu Mar 26 20:54:54 2015
 * 
 * @brief    compressed sparse row data structure a good data structure for a fixed topology 
 * 
 * 
 */

#ifndef COMPOSLIST_H
#define COMPOSLIST_H

#include<cstdio>
#include<vector>
#include<cassert>
#include<algorithm>
#include"heap.h"
#include<queue>
#include<set>
#include<limits>
#include<iostream>
using std::vector;
using namespace std;



template <class E, class W, class C>
struct endElement{
  C captical;
  W weight;
  E snk;
  endElement(  ):captical( 0 ),  weight( 0 ), snk( 0 ){
    
  }
};

template<class E>
struct startElement{
  E link;
  E src;
  startElement(  ):link( 0 ), src( 0 ){
    
  }
  
};

template<class E, class W>
struct precedence
{
  E link;
  W weight;
  E vertex;
  precedence(  ):link( 0 ),  vertex( 0 ){
  }
  precedence & operator=( const precedence&  other){
    link=other.link ;
    weight=other.weight;
    vertex=other.vertex;
    return *this;
  }

};

template <typename E>
struct tempElment
{
  E id;
  E src;
};

template<typename E>
static bool tempElmentCmp(const tempElment<E> &lhs, const tempElment<E> &rhs  ) {
  return lhs.src<rhs.src;
}



template<class T>
struct LESSOR{
  bool operator(  )(const T &x, const T &y ) const{
    return x.first< y.first;
  }
  
};


template<typename E=int, typename W=float, typename C=float>
class compressed_sparse_row_graph{

private:
  E vertex_num;
  E link_num;
  //out links
  vector<E> outIndex;
  vector<endElement<E,W,C> > link_ends;
  //in links
  vector<E> inIndex;
  vector<startElement<E> > link_starts;

  vector<E> link_map;// link_ends  --> orignal link
  vector<E> reLink_map; // orginal link --> link_end
  vector<precedence<E,W> > shortPaths;  // allpair shortest path

  W infi_value;
  int _findSrc( const E link, E & src )const {

    size_t  start, end, mid;

    src=vertex_num+1;
    if(link >=link_num ) return 0;

    start=0;
    end=vertex_num-1;
    while( end> start ){

      mid=( start+end )/2;
      if(link< outIndex[ mid ] ){
        end=mid-1;
      }
      else if( link>= outIndex[ mid+1 ] ){
        start=mid+1;
      }
      else  {
        src=mid;
        return 1;
      }
    }


    if( link>= outIndex[ start ] && link<  outIndex[ start+1 ] ){
      src=start;
      return 1;    
    }

    src=end;
    return 1;
    
  }
  int _findSnk( const E link, E & snk ) const{
    snk=vertex_num+1;
    if( link>=link_num ) return 0;
    snk=link_ends[ link ].snk;
    return 1;

  }

  inline int _getSrcSnk( const E link, E &src, E &snk  ) const{
    if(  0== _findSrc( link, src )) return 0;
    _findSnk( link, snk );
    return 1;
  }

  int _getShortPath( const E src, const E snk, vector<E> & path ) const{
    assert( src< vertex_num && snk < vertex_num );
    assert( src!=snk );

    const size_t shift= src*vertex_num;
    if(shortPaths[ shift +snk ].vertex> vertex_num )
      return  0;
  
    E  current =snk;
    while( current!=src ){
      path.push_back( current );
      path.push_back( shortPaths[ shift +current ].link   );
      current=shortPaths[ shift +current ].vertex;
    }
    path.push_back( current );
    std::reverse(path.begin(  ), path.end(  ) );
    return 1;
    

  }
  int _increaseLinkWeight(const  E link, const W inc ){
    
    assert( link< link_num );
    assert( inc>=0 );
    if( 0==inc ) return 0;

    link_ends[ link ].weight+=inc;

    size_t i;
    size_t src, snk;
    double orignalW, weight;
    _getSrcSnk( link, src, snk );

#pragma omp parallel for 
    for (i = 0; i < vertex_num; i++) {
      size_t j, current, inDegree, outDegree;
      const  size_t shift=i*vertex_num;
      if(link==shortPaths[ shift+snk ].link  ){
        queue<size_t>  wait;
        set<size_t> pass;
        pass.insert( i );
        wait.push(snk );
        pass.insert( snk );
        while (!wait.empty(  )) {
          current=wait.front(  );
          wait.pop(  );
          orignalW=shortPaths[shift+snk].weight;
          shortPaths[shift+snk].weight=infi_value;
          precedence<E,W> &dummy_current=shortPaths[shift + current];

          inDegree=getInDegree( current );
          for( j=0; j< inDegree; j++ ){
            const startElement<E> &neighbour=link_starts[ inIndex[ current ]+j ];

            weight=shortPaths[ shift+ neighbour.src ].weight+ link_ends[ neighbour.link ].weight;
            if( weight< dummy_current.weight ){
              dummy_current.link=  neighbour.link;
              dummy_current.weight=weight;
              dummy_current.vertex=neighbour.src;
            
            }
          }
          if( dummy_current.weight> orignalW ){
            outDegree=getOutDegree( current );
            for( j=0; j< outDegree; j++ ){
              if(  pass.find(link_ends[outIndex[ current ]+j].snk)== pass.end(  )  ){
                pass.insert( link_ends[outIndex[ current ]+j ].snk );
                wait.push( link_ends[outIndex[ current ]+j ].snk );
              
              }
            
            }

          }
        
        }

      }
    
    }
    return 1;
  }


public:
  
  compressed_sparse_row_graph( ):vertex_num( 0 ), link_num( 0 ){
  }

  void initial(vector<E> &srcs, vector< E > &snks, vector<W> &weights  ){

    assert(srcs.size( )==  snks.size(  ) && srcs.size(  )== weights.size(  ) );
    if ( 0== srcs.size(  ) ) return ;
    link_num=srcs.size(  );
    link_map.resize(link_num  );
    reLink_map.resize( link_num );

    vector<tempElment<E> > tContian;
   

    E i, j, ver_max;
    tempElment<E> dummy_pred;

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
    precedence<E,W> inifPre;
    inifPre.link=link_num+1;
    inifPre.weight=numeric_limits<W>::max( )/10;
    inifPre.vertex=vertex_num+10;      
    shortPaths.resize(  vertex_num*vertex_num, inifPre);

    for(i=0 ; i< srcs.size( ) ; i++ ){
      dummy_pred.id=i;
      dummy_pred.src=srcs[ i ];
      tContian.push_back( dummy_pred );
    }
    std::sort( tContian.begin(  ), tContian.end(  ), tempElmentCmp<E>  );

    outIndex.push_back( 0 );
    i=0;
    while(i< tContian[ 0 ].src  ){
      outIndex.push_back( 0 );
      i++;
    }

    endElement<E,W,C> temp;
    temp.weight=weights[tContian[ 0 ].id ];
    temp.snk=snks[tContian[ 0 ].id ];


    link_ends.push_back( temp );
    link_map[ 0 ]=tContian[ 0 ].id;
    reLink_map[ tContian[ 0 ].id ]=0;
  
    for( i=1; i< tContian.size( ) ;i++ ){
      if(  tContian[ i ].src!= tContian[ i-1 ].src ){
        for(j=tContian[ i-1 ].src; j< tContian[ i ].src; j++ ){
          outIndex.push_back(link_ends.size( ) );
        }
      }
      temp.weight=weights[ tContian[ i ].id ];
      temp.snk=snks[ tContian[ i ].id ];


      link_ends.push_back( temp );
      link_map[ i ]=tContian[ i ].id;
      reLink_map[tContian[ i ].id  ]=i;
    }

    for( j=tContian[ i-1 ].src; j< vertex_num; j++ ){
      outIndex.push_back(link_ends.size( ) );
    }


  
    tContian.clear(  );
    for(i=0 ; i< srcs.size( ) ; i++ ){
      dummy_pred.id=i;
      dummy_pred.src=snks[ i ];
      tContian.push_back( dummy_pred );
    }

    std::sort( tContian.begin(  ), tContian.end(  ), tempElmentCmp<E>  );

    inIndex.push_back( 0 );
    i=0;
    while(i< tContian[ 0 ].src  ){
      inIndex.push_back( 0 );
      i++;
    }
    startElement<E> dummy;
    dummy.link=tContian[ 0 ].id;
    dummy.src=srcs[ tContian[ 0 ].id ];

    link_starts.push_back( dummy );
  
    for (i = 1; i < tContian.size(  ); i++) {

      if( tContian[ i ].src!= tContian[ i-1 ].src ){
        for( j=tContian[ i-1 ].src; j< tContian[ i ].src; j++ ){
          inIndex.push_back( link_starts.size(  ) );        
        }
      }
      dummy.link=tContian[ i ].id;
      dummy.src=srcs[ tContian[ i ].id ];
      link_starts.push_back( dummy );
    }

    for( j=tContian[ i-1 ].src; j< vertex_num; j++ ){
      inIndex.push_back(link_starts.size( ) );
    }


    compute_allPair_shortest_path(  );
    

  }
  inline size_t getVertex_num( void  )const{
    return vertex_num;
  }

  void setInfi( const E infi ){
    infi_value=infi;
  }

  inline size_t getLink_num ( void  ) const{
    return link_num;
  }

  inline int getOutDegree(E vertex) const{
    assert( vertex< vertex_num );
    return outIndex[ vertex+1 ]-outIndex[ vertex ];
  }
  inline int getInDegree( E vertex) const {
    assert( vertex< vertex_num );
    return inIndex[vertex+1 ]-inIndex[ vertex ];

  }
  inline endElement<E,W,C> & getLink(E vertex, int k ) {
    assert(  vertex< vertex_num && k< getOutDegree( vertex )  );
    return  link_ends[outIndex[vertex ]+k ];
  }
  inline W getWeight( const E link ) const{
    return link_ends[reLink_map[ link ]    ].weight;
  }
  inline C getCapitial(const E link ) const{
    return link_ends[reLink_map[ link ]    ].captical;
  }

  inline int findSrc( const E link, E & src )const {
    return _findSrc( reLink_map[ link ], src );
  }
  inline int findSnk( const E link, E & snk ) const{
    return _findSnk(reLink_map[ link ], snk  );
  }

  inline int getSrcSnk( const E link, E &src, E &snk  ) const{
    return _getSrcSnk(reLink_map[link], src, snk  );
  }

  inline int increaseLinkWeight(const  size_t link, const W inc ){
    return _increaseLinkWeight(reLink_map[link], inc ) ;

  }
  inline void setLinkWeight( const E link, const E weight ){
    link_ends[reLink_map[ link ]  ].weight=weight;
  }

  inline void setLinkCap( const E link, const C captical  ){
    link_ends[reLink_map[ link ]  ].captical=captical;
  }

  
  void  compute_sourceallPair_shortest_path_dijkstra(const E src ,   bool resize=true ){

    typedef pair<W ,E> PII;
    const  size_t shift=src*vertex_num;
    size_t j, current, outDegree ;
    W weight;
    //    Q.clear(  );
    if( resize ){

      precedence<E,W> inifPre;
      inifPre.link=link_num+1;
      inifPre.weight=numeric_limits<W>::max( )/10;
      inifPre.vertex=vertex_num+10;  
      fill( shortPaths.begin(  ), shortPaths.begin(  )+vertex_num , inifPre);

    }

    shortPaths[shift+src].weight=0;
    //    vector<precedence<E,W> > tempshortPaths(shortPaths);


    LESSOR<PII> order;
    Fixed_heap<W,E , LESSOR<PII> > Q( order, vertex_num );
   //    Q.setCap(vertex_num  );
    //priority_queue<PII, vector<PII>, greater<PII> >Q;
    
    Q.push( make_pair( 0.0, src ) );
    while( !Q.empty(  ) ){
      PII p=Q.top(  );
      Q.pop(  );
      current=p.second;
      outDegree=getOutDegree( current );

      for(j=0; j< outDegree; j++  ){
        const endElement<E,W,C> &neighbour=link_ends[outIndex[current]+j];
        precedence<E,W> &dummy_pred=shortPaths[shift + neighbour.snk];
        weight=shortPaths[ shift+current ].weight+ neighbour.weight;
        if( weight <dummy_pred.weight   ){
          dummy_pred.link= outIndex[current]+j;
          dummy_pred.weight =weight;
          dummy_pred.vertex=current;
          Q.push(make_pair(weight,  neighbour.snk) );


        }

      }
    }    


    // priority_queue<PII, vector<PII>, greater<PII> >Q1;
    // Q1.push( make_pair( 0.0, src ) );

    // while( !Q1.empty(  ) ){
    //   PII p=Q1.top(  );
    //   Q1.pop(  );
    //   current=p.second;
    //   outDegree=getOutDegree( current );

    //   for(j=0; j< outDegree; j++  ){
    //     const endElement<E,W,C> &neighbour=link_ends[outIndex[current]+j];
    //     precedence<E,W> &dummy_pred=tempshortPaths[shift + neighbour.snk];
    //     weight=tempshortPaths[ shift+current ].weight+ neighbour.weight;
    //     if( weight <dummy_pred.weight   ){
    //       dummy_pred.link= outIndex[current]+j;
    //       dummy_pred.weight =weight;
    //       dummy_pred.vertex=current;
    //       Q1.push(make_pair(weight,  neighbour.snk) );


    //     }

    //   }
    // }    

    // for (j = 0; j < vertex_num; j++) {
      
    //   if( shortPaths[ shift+j ].link!= tempshortPaths[ shift+j ].link ||  shortPaths[ shift+j ].vertex!= tempshortPaths[ shift+j ].vertex)  
    //     std::cout << "error" << std::endl;
  
    // }

  }
  /** 
   * http://en.wikipedia.org/wiki/Suurballe%27s_algorithm
   * 
   * @param src 
   * @param snk 
   * @param path1 
   * @param path2 
   */
  void suurballe_shortest(  const E src, const E snk,  vector<E> & path1,  vector<E> & path2)
  {
    const  size_t shift=src*vertex_num;

    /**
     *  fisrt dijkstra shortest path
     * 
     */

    compute_sourceallPair_shortest_path_dijkstra(src  );

    vector<E> tempp;
    if(_getShortPath( src, snk, tempp ) <=0 ){
      return ;
    }
    vector<W>  newW(link_num  );
    E i;
    E srcc, snkk;
    /**
     *  new edge weight
     * 
     */

    for( i=0; i< link_num; i++ )
      {
        _getSrcSnk(  i, srcc, snkk);
        newW[ i ]=link_ends[ i ].weight+shortPaths[ shift+srcc ].weight-shortPaths[ shift+snkk ].weight;
        assert( newW[ i ] )>=0);
  }

  vector<E> path1link;
  i=0;
  while (i< tempp.size(  ) -1) {
    path1link.push_back(tempp[ i+1 ]  );
    i+=2;
  }

  /**
   *  second dijkstra shortest path
   * 
   */

  typedef pair<W ,E> PII;

  E link;
  size_t j, current, outDegree ;
  W weight;

  vector<W> dis(vertex_num, numeric_limits<W>::max( )/10 );
  vector<E> parent(vertex_num, vertex_num+10  );
  vector<E> shortp(vertex_num, link_num+10  );
    
  dis[ src ]=0;

  LESSOR<PII> order;
  Fixed_heap<W,E , LESSOR<PII> > Q( order, vertex_num );

    
  Q.push( make_pair( 0.0, src ) );
  while( !Q.empty(  ) ){
    PII p=Q.top(  );
    current=p.second;
    if( snk==current ){
      break;
    }
    Q.pop(  );
      
    outDegree=getOutDegree( current );

    for(j=0; j< outDegree; j++  ){
      link=outIndex[current]+j;
      // reverse all directtion of link in path1
      if(find( path1link.begin(  ), path1link.end(  ), link )!= path1link.end(  )  )
        continue;

      const endElement<E,W,C> &neighbour=link_ends[outIndex[current]+j];

      if( src==neighbour.snk ) continue; // cut all link direct to src

      weight=dis[ current ] + newW[ link ];

      if( weight <dis[neighbour.snk  ]   ){
        shortp[ neighbour.snk ]= link;
        dis[ neighbour.snk ]=weight;
        parent[ neighbour.snk ]=current;

        Q.push(make_pair(weight,  neighbour.snk) );

      }
    }
      
    // reverse all directtion of link in path1 and  cut all link direct to src
    for(j=0; j< path1link.size(  );  j++ )
      {
        link=path1link[ j ];
        _getSrcSnk( link, srcc, snkk  ); 
        if( current==snkk && srcc!= src)
          {
            weight=dis[ current ]+newW[ link ];
            assert( 0.0== newW[ link ]);
            if( weight< dis[ srcc ] ){
              shortp[ srcc ]=link;
              dis[ srcc ]=weight;
              parent[ srcc ]= current;
              Q.push( make_pair( weight, srcc ) );

            }
          }
          
      }

  
    if( Q.empty(  ) ) 
      return ;

    vector<E>  path2link;
    current=snk;
    while (snk!=src) {
      path2link.push_back( shortp[ current ] );
      current=parent[current  ];
    }

    /**
     * delete the same links in path1link and path2link
     * 
     */
    
    sort(path2link.begin(  ) , path2link.end(  )  );
    sort( path1link.begin(  ), path1link.end(  ) );
    vector<E> sameLink( path1link.size(  ) );
    vector<E>::iterator it=set_intersection( path1link.begin(  ), path1link.end(  ), path2link.begin(  ), path2link.end(  ), sameLink.begin(  ) );
    sameLink.resize( it-sameLink.begin(  ) );
    
    for( j=0; j< sameLink.size(  ); j++ ){
      it=find(path1link.begin(  ), path1link.end(  ), sameLink[ j ]);
      if( it!=path1link.end(  ) ){
        path1link.erase( it );
      }

    }
    
    for( j=0; j< sameLink.size(  ); j++ ){
      it=find(path2link.begin(  ), path2link.end(  ), sameLink[ j ]);
      if( it!=path2link.end(  ) ){
        path2link.erase( it );
      }
    }

    path1link.insert( path1link.end(  ),  path2link.begin(  ), path2link.end(  ) );



    /**
     * obtain the two disjoint paths
     * 
     */

    current=src;
    while (current!=snk) {
      j=0;    
      for( j=0;  j< path1link.size(  ) ; j++)
        {
          link=path1link[ j ];
          _findSrc( path1link[ j ] , srcc);
          if( srcc==current ) break;
        }

      assert( j<path1link.size(  )  );

      path1link.erase( path1link.begin(  )+j );

      path1.push_back( link );
      _findSnk( link, snkk );
      current=snkk;

    }


    current=src;
    while (current!=snk) {
      j=0;    
      for( j=0;  j< path1link.size(  ) ; j++)     {
        link=path1link[ j ];
        _findSrc( path1link[ j ] , srcc);
        if( srcc==current ) break;
      }

      assert( j<path1link.size(  )  );

      path1link.erase( path1link.begin(  )+j );

      path2.push_back( link );
      _findSnk( link, snkk );
      current=snkk;

    }
    
    
  }



  void compute_sourceallPair_shortest_path(const E src  ){

    const  size_t shift=src*vertex_num;
    size_t j, current, outDegree ;
    W weight;
    shortPaths[shift+src].weight=0;
    queue<E>  wait;
    wait.push(src );
    while( !wait.empty(  ) ){
      current=wait.front(  );
      wait.pop(  );
      outDegree=getOutDegree( current );
      for( j=0; j< outDegree; j++ ){
        const    endElement<E,W,C> &neighbour=link_ends[outIndex[current]+j];
        precedence<E,W> &dummy_pred=shortPaths[shift + neighbour.snk];
        weight=shortPaths[ shift+current ].weight+ neighbour.weight;

        if( weight <dummy_pred.weight   ){
          dummy_pred.link= outIndex[current]+j;
          dummy_pred.weight =weight;
          dummy_pred.vertex=current;
          wait.push( neighbour.snk );
        }
      }
    }  
    

  }

  void compute_allPair_shortest_path(  ){
    
    E i;
    typedef pair<W ,E> PII;
    precedence<E,W> inifPre;
    inifPre.link=link_num+1;
    inifPre.weight=numeric_limits<W>::max( )/10;
    inifPre.vertex=vertex_num+10;

    //    shortPaths.resize(vertex_num*vertex_num);
    fill( shortPaths.begin(  ), shortPaths.end(  ) , inifPre);
    // LESSOR<PII> order;
    // Fixed_heap<W,E , LESSOR<PII> > Q( order, vertex_num );

#pragma omp parallel for
    for( i=0; i< vertex_num; i++ ){
      compute_sourceallPair_shortest_path_dijkstra( i ,    false);
    }
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
  
  int  getShortPath( const E src, const E snk, vector<E> & path ) const {
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
 

  static void printPath( vector<E>&path ){
    if(0== path.size(  ) ) return;
    
    size_t i=0;
    string ss;
    
    
    

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


