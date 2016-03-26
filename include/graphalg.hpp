/**
 * @file   graphalg.hpp
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Sat Mar 26 23:42:07 2016
 *
 * @brief  graph common algorithm
 *
 *
 */

#ifndef __GRAPH_ALG_H
#define __GRAPH_ALG_H

namespace fast_graph {
template <typename T>
struct LESSOR_T {
  bool operator()(const T &x, const T &y) const { return x.first < y.first; }
};

template <typename WV,  typename W>
W  path_cost( const WV& NW, const vector<int> &path, W  ){
  W re=0;
  for( vector<int>::const_iterator it=path.begin(  ); it!= path.end( ); it++ ){
    re+=NW[ *it ];
  }
  return re;
}


template <typename G, typename WV,  typename W>
bool dijkstra_shortest_path(const G &g, const WV &NW, const W inif,
                            const int src, const int snk, vector<int> &path) {
  typedef pair<W, int> PII;
  path.clear();
  if (src < 0 || snk < 0) return false;
  if (src == snk) return true;
  int vertex_num = g.getVertex_num();
  if (src >= vertex_num || snk >= vertex_num) return false;

  vector<int> preLink(vertex_num, -1);
  vector<W> dis(vertex_num, inif  );
  LESSOR_T<PII> order;
  size_t j, outDegree;
  int link, next;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII> > Q(order, vertex_num);
  dis[ src ]=0;
  Q.push(make_pair(0.0, src));
  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;
    if (current == snk) {
      while (current != src) {
        path.push_back(preLink[current]);
        g.findRhs(preLink[current], current, current  );

      }
      reverse(path.begin(), path.end());
      return true;
    }
    Q.pop();


    outDegree = g.getOutDegree(current);
    W current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = g.getAdj( current, j );

      weight = current_weight + NW[ link ];
      
      g.findRhs( link, current, next );

      if (weight < dis[snk] && weight < dis[next]) {
        dis[ next ]=weight;
        preLink[next] = link;
        Q.push(make_pair(weight, next));
      }
    }
  }
  return false;
}
}

#endif
