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

#include <algorithm>
#include <limits>
#include<map>
#include <queue>
#include <stack>
#include <unordered_set>
#include <utility>
#include <vector>
#include "graph.h"
#include "heap.hpp"
#include <mutex> 

#if  (defined(__APPLE__) && defined(__MACH__))

#else
#include <omp.h>

#endif

namespace raptor {
using std::vector;
using std::unordered_set;

using std::queue;
using std::stack;
using std::pair;
using namespace std;


template <typename T>
struct LESSOR_T {
  bool operator()(const T &x, const T &y) const { return x.first < y.first; }
};
template<typename W>
struct Data{
  int vertex_num;
  typedef pair<W, int> PII;
  vector<int> preLink;
  vector<W> dis;
  vector<char> check;
  Fixed_heap<W, int, LESSOR_T<PII>> Q;
  vector<int> bpreLink;
  vector<W> bdis;
  Fixed_heap<W, int, LESSOR_T<PII>> bQ;
  LESSOR_T<PII> order;

  Data():vertex_num(0),Q(order), bQ(order){
  }
  Data(int num, W inf):vertex_num(num),preLink(num, -1), dis(num, inf),
                         check( num, 0), Q(order, num), bpreLink(num, -1),
                         bdis(num, inf), bQ(order, num) {
  }
  void reset(W inf) {
    int temp;
    const vector<int>& qpassNodes=Q.getPassNodes();
    assert(qpassNodes.size()<= preLink.size());
    for(vector<int>::const_iterator it=qpassNodes.begin(); it!=qpassNodes.end(); it++){
      temp=*it;
      preLink[temp]=-1;
      check[temp]=0;      
    }
    const vector<int>& bqpassNodes=bQ.getPassNodes();
    for(vector<int>::const_iterator it= bqpassNodes.begin(); it!= bqpassNodes.end(); it++){
      temp=*it;
      bpreLink[temp]=-1;
      check[temp]=0;
    }
    fill(dis.begin(), dis.end(), inf);
    fill(bdis.begin(), bdis.end(), inf);
    Q.clear();
    bQ.clear();
  }
  
};
/** 
 * 
 * @brief has some when  parallel invoke
 * @param vertex_num 
 * @param inf 
 * @param id 
 * @param isGet 
 * 
 * @return 
 */

template<typename W>
static Data<W>& getData(int vertex_num, W inf, int &id, bool isGet=true){
  typedef pair<W, int> PII;
  
  static  std::vector< Data<W> > datas;
 
  static  vector<char> isUsed;
  

  static std::mutex mtx;



  #pragma omp critical
  {
    if(isGet)
    {
      id=-1;
      for (size_t i=0; i< datas.size(); i++){
        if(!isUsed[i] &&  datas[i].vertex_num>=vertex_num){
          id=i;
          isUsed[id]=1;
          datas[id].reset(inf);
          break;
        }
      }
      if(id<0){
        Data<W>  d(vertex_num,  inf );
        id=datas.size();
        int thread_num=8;
        for(int i=0; i< thread_num; i++){
          datas.push_back(d);
          isUsed.push_back(0);          
        }
      }
    }else{
      isUsed[id]=0;
    }
  }
  
  return datas[id];  
}

template <typename WV, typename W>
W path_cost(const WV &NW, const vector<int> &path, W) {
  W re = 0;
  for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
    re += NW[*it];
  }
  return re;
}
template <typename G>
bool isValidatePath(const G &graph, const int &src, const int &snk,
                    const vector<int> &path) {
  int current = src;
  for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
    if (!graph.findRhs(*it, current, current)) {
      return false;
    }
  }
  return current == snk;
}

template <typename G>
bool isSimplePath(const G &graph, const int src, const int snk,
                  const vector<int> &path) {
  int current = src;
  vector<int> nodes;
  nodes.push_back(current);
  for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
    if (!graph.findRhs(*it, current, current)) {
      return false;
    }
    nodes.push_back(current);
  }
  if (current != snk) {
    return false;
  }

  sort(nodes.begin(), nodes.end());
  int last = -1;
  for (vector<int>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    if (*it == last) {
      return false;
    }
    last = *it;
  }
  return true;
}

/**
 * @brief https://en.wikipedia.org/wiki/DFS
 *
 * @param graph
 * @param src
 * @param pass_nodes
 *
 * @return
 */
template <typename G>
void dfs_search(const G &graph, const int src, vector<int> &pass_nodes) {
  pass_nodes.clear();
  int vertex_num = graph.getVertex_num();
  if (src < 0) {
  }

  if (src >= vertex_num) {
  }

  int current, j, link, outDegree, tempSnk;

  vector<int> child_id(vertex_num, -1);
  pass_nodes.push_back(src);
  child_id[src] = 0;
  stack<int> Q;
  Q.push(src);

  while (!Q.empty()) {
    current = Q.top();
    outDegree = graph.getOutDegree(current);
    if (outDegree == child_id[current]) {
      Q.pop();

    } else {
      int k = child_id[current];
      child_id[current]++;
      link = graph.getAdj(current, k);
      graph.findRhs(link, current, tempSnk);
      if (child_id[tempSnk] < 0) {
        child_id[tempSnk] = 0;
        Q.push(tempSnk);
        pass_nodes.push_back(tempSnk);
      }
    }
  }
}

/**
 * @brief https://en.wikipedia.org/wiki/Breadth-first_search
 *
 * @param graph
 * @param src
 * @param pass_nodes
 */
template <typename G>
void bfs_search(const G &graph, const int src, vector<int> &pass_nodes) {
  pass_nodes.clear();
  int vertex_num = graph.getVertex_num();
  if (src < 0) {
  }

  if (src >= vertex_num) {
  }

  int current, j, link, outDegree, tempSnk;
  vector<int> Q;
  vector<bool> pass(vertex_num, false);
  pass[src] = true;
  while (!Q.empty()) {
    vector<int> secondQ;
    for (vector<int>::iterator it = Q.begin(); it != Q.end(); it++) {
      current = *it;
      outDegree = graph.getOutDegree(current);
      for (j = 0; j < outDegree; j++) {
        link = graph.getAdj(current, j);
        graph.findRhs(link, current, tempSnk);
        if (!pass[tempSnk]) {
          pass[tempSnk] = true;
          secondQ.push_back(tempSnk);
        }
      }
    }
    Q = secondQ;
  }
}

template <typename G, typename WV, typename W>
bool dijkstra_shortest_path(const G &graph, const WV &NW, const int src,
                            const int snk, vector<int> &path, const W inf) {
  typedef pair<W, int> PII;
  path.clear();
  if (src < 0 || snk < 0) {
    return false;
  }
  if (src == snk) {
    return true;
  }
  int vertex_num = graph.getVertex_num();
  if (src >= vertex_num || snk >= vertex_num) {
    return false;
  }

  vector<int> preLink(vertex_num, -1);
  vector<W> dis(vertex_num, inf);
  LESSOR_T<PII> order;
  size_t j, outDegree;
  int link, tempSnk;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
  dis[src] = 0;
  Q.push(make_pair(0.0, src));
  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;
    if (current == snk) {
      while (current != src) {
        path.push_back(preLink[current]);
        graph.findRhs(preLink[current], current, current);
      }
      reverse(path.begin(), path.end());
      return true;
    }
    Q.pop();

    outDegree = graph.getOutDegree(current);
    W current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = graph.getAdj(current, j);

      weight = current_weight + NW[link];

      graph.findRhs(link, current, tempSnk);

      if (weight < dis[snk] && weight < dis[tempSnk]) {
        dis[tempSnk] = weight;
        preLink[tempSnk] = link;
        Q.push(make_pair(weight, tempSnk));
      }
    }
  }
  return false;
}

/**
 * @brief https://en.wikipedia.org/wiki/A*_search_algorithm
 *
 * @param graph
 * @param NW
 * @param h  heuristic function
 * @param src
 * @param snk
 * @param path
 * @param inf
 *
 * @return
 */

template <typename G, typename WV, typename W, typename H>
bool astar_shortest_path(const G &graph, const WV &NW, const H &h,
                         const int src, const int snk, vector<int> &path,
                         const W inf) {
  typedef pair<W, int> PII;
  path.clear();
  if (src < 0 || snk < 0) {
    return false;
  }
  if (src == snk) {
    return true;
  }
  int vertex_num = graph.getVertex_num();
  if (src >= vertex_num || snk >= vertex_num) {
    return false;
  }

  vector<int> preLink(vertex_num, -1);
  vector<W> dis(vertex_num, inf);
  LESSOR_T<PII> order;
  size_t j, outDegree;
  int link, tempSnk = 0;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);

  dis[src] = 0;
  Q.push(make_pair(0.0, src));
  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;
    if (current == snk) {
      while (current != src) {
        path.push_back(preLink[current]);
        graph.findRhs(preLink[current], current, current);
      }
      reverse(path.begin(), path.end());
      return true;
    }
    Q.pop();

    outDegree = graph.getOutDegree(current);
    W current_weight = dis[current];
    for (j = 0; j < outDegree; j++) {
      link = graph.getAdj(current, j);

      weight = current_weight + NW[link];

      graph.findRhs(link, current, tempSnk);

      if (weight < dis[snk] && weight < dis[tempSnk]) {
        dis[tempSnk] = weight;
        preLink[tempSnk] = link;
        Q.push(make_pair(weight + h[tempSnk], tempSnk));
      }
    }
  }
  return false;
}

template <typename G, typename WV, typename W>
void dijkstra_shortest_tree(const G &graph, const WV &NW, const int src,
                            vector<int> &preLink, const W inf) {
  typedef pair<W, int> PII;
  int vertex_num = graph.getVertex_num();
  preLink.resize(vertex_num);
  fill(preLink.begin(), preLink.end(), -1);

  if (src < 0) {
    return;
  }

  if (src >= vertex_num) {
    return;
  }

  vector<W> dis(vertex_num, inf);
  LESSOR_T<PII> order;
  size_t j, outDegree;
  int link, tempSnk;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
  dis[src] = 0;
  Q.push(make_pair(0.0, src));
  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;

    Q.pop();

    outDegree = graph.getOutDegree(current);
    W current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = graph.getAdj(current, j);

      weight = current_weight + NW[link];

      graph.findRhs(link, current, tempSnk);

      if (weight < dis[tempSnk]) {
        dis[tempSnk] = weight;
        preLink[tempSnk] = link;
        Q.push(make_pair(weight, tempSnk));
      }
    }
  }
}

template <typename G, typename WV, typename W>
void dijkstra_shortest_tree(const G &graph, const WV &NW, const int src,
                            vector<int> &preLink, vector<W> &dis, const W inf) {
  typedef pair<W, int> PII;
  int vertex_num = graph.getVertex_num();
  preLink.resize(vertex_num);
  dis.resize(vertex_num);
  fill(preLink.begin(), preLink.end(), -1);
  fill(dis.begin(), dis.end(), inf);
  if (src < 0) {
    return;
  }

  if (src >= vertex_num) {
    return;
  }

  LESSOR_T<PII> order;
  size_t j, outDegree;
  int link, tempSnk;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
  dis[src] = 0;
  Q.push(make_pair(0.0, src));
  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;

    Q.pop();

    outDegree = graph.getOutDegree(current);
    W current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = graph.getAdj(current, j);

      weight = current_weight + NW[link];

      graph.findRhs(link, current, tempSnk);

      if (weight < dis[tempSnk]) {
        dis[tempSnk] = weight;
        preLink[tempSnk] = link;
        Q.push(make_pair(weight, tempSnk));
      }
    }
  }
}

template <typename G, typename WV, typename W>
void dijkstra_shortest_retree(const G &graph, const WV &NW, const int snk,
                              vector<int> &preLink, vector<W> &dis,
                              const W inf) {
  typedef pair<W, int> PII;
  int vertex_num = graph.getVertex_num();
  preLink.resize(vertex_num);
  dis.resize(vertex_num);
  fill(preLink.begin(), preLink.end(), -1);
  fill(dis.begin(), dis.end(), inf);
  if (snk < 0) {
    return;
  }

  if (snk >= vertex_num) {
    return;
  }

  LESSOR_T<PII> order;
  size_t j, inDegree;
  int link, tempSnk = snk;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
  dis[snk] = 0;
  Q.push(make_pair(0.0, snk));
  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;

    Q.pop();

    inDegree = graph.getInDegree(current);
    W current_weight = p.first;
    for (j = 0; j < inDegree; j++) {
      link = graph.getReAdj(current, j);

      weight = current_weight + NW[link];

      graph.findRhs(link, current, tempSnk);

      if (weight < dis[tempSnk]) {
        dis[tempSnk] = weight;
        preLink[tempSnk] = link;
        Q.push(make_pair(weight, tempSnk));
      }
    }
  }
}

template <typename G, typename WV, typename W>
void dijkstra_shortest_retree(const G &graph, const WV &NW, const int snk,
                              vector<W> &dis, const W inf) {
  vector<int> preLink;
  dijkstra_shortest_retree(graph, NW, snk, preLink, dis, inf);
}

template <typename G, typename WV, typename W>
bool bidijkstra_shortest_path(const G &graph, const WV &NW, const int src,
                              const int snk, vector<int> &path, const W inf) {
  typedef pair<W, int> PII;
  path.clear();
  if (src < 0 || snk < 0) {
    return false;
  }
  if (src == snk) {
    return true;
  }
  int vertex_num = graph.getVertex_num();
  if (src >= vertex_num || snk >= vertex_num) {
    return false;
  }
#ifdef STATIC_TABLE
  int did;
  Data<W> &data=getData(vertex_num, inf, did);
  
  vector<int>& preLink=data.preLink;
  vector<W>& dis=data.dis;
  vector<char>& check=data.check;
#else


  vector<int> preLink(vertex_num, -1);
  vector<W> dis(vertex_num,inf);
  vector<char> check(vertex_num, 0);
#endif

  LESSOR_T<PII> order;
  size_t j, outDegree, inDegree;
  int link, tempSnk = 0;
  int current;
  W weight;
#ifdef STATIC_TABLE
  Fixed_heap<W, int, LESSOR_T<PII>>& Q=data.Q;
#else
  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
#endif

  dis[src] = 0;
  preLink[src]=0;
  Q.push(make_pair(0.0, src));
#ifdef STATIC_TABLE
  vector<int>& bpreLink=data.bpreLink;
  vector<W>& bdis=data.bdis;
#else 
  vector<int> bpreLink(vertex_num, -1);
  vector<W> bdis(vertex_num, inf);
#endif
  

#ifdef STATIC_TABLE
  Fixed_heap<W, int, LESSOR_T<PII>>& bQ=data.bQ;
#else 
  Fixed_heap<W, int, LESSOR_T<PII>> bQ(order, vertex_num);
#endif
  bdis[snk] = 0;
  bpreLink[snk]=0;
  bQ.push(make_pair(0.0, snk));
  bool re = false;
  W best_dis = inf;
  int best_current = -1;
  if (graph.isDirect()) {
    while (!Q.empty() && !bQ.empty()) {
      PII p = Q.top();
      current = p.second;
      if (check[current]) {
        re = true;
        break;
      }
      check[current] = 1;
      Q.pop();

      outDegree = graph.getOutDegree(current);
      W current_weight = p.first;
      for (j = 0; j < outDegree; j++) {
        link = graph.getAdj(current, j);
        weight = current_weight + NW[link];
        graph.findSnk(link, tempSnk);
        if (weight < dis[tempSnk]) {
          dis[tempSnk] = weight;
          preLink[tempSnk] = link;
          Q.push(make_pair(weight, tempSnk));
        }
      }
      p = bQ.top();
      current = p.second;
      if (check[current]) {
        re = true;
        break;
      }
      check[current] = 1;
      bQ.pop();
      inDegree = graph.getInDegree(current);
      current_weight = bdis[current];
      for (j = 0; j < inDegree; j++) {
        link = graph.getReAdj(current, j);
        weight = current_weight + NW[link];
        graph.findSrc(link, tempSnk);
        if (weight < bdis[tempSnk]) {
          bdis[tempSnk] = weight;
          bpreLink[tempSnk] = link;
          bQ.push(make_pair(weight, tempSnk));
        }
      }
    }
  } else {
    while (!Q.empty() && !bQ.empty()) {
      PII p = Q.top();
      current = p.second;
      if (check[current]) {
        re = true;
        break;
      }
      check[current] = 1;
      Q.pop();

      outDegree = graph.getOutDegree(current);
      W current_weight = p.first;
      for (j = 0; j < outDegree; j++) {
        link = graph.getAdj(current, j);

        weight = current_weight + NW[link];

        graph.findRhs(link, current, tempSnk);

        if (weight < dis[tempSnk]) {
          dis[tempSnk] = weight;
          preLink[tempSnk] = link;
          Q.push(make_pair(weight, tempSnk));
        }
      }
      p = bQ.top();
      current = p.second;
      if (check[current]) {
        re = true;
        break;
      }
      check[current] = 1;
      bQ.pop();
      inDegree = graph.getInDegree(current);
      current_weight = bdis[current];
      for (j = 0; j < inDegree; j++) {
        link = graph.getReAdj(current, j);
        weight = current_weight + NW[link];
        graph.findRhs(link, current, tempSnk);
        if (weight < bdis[tempSnk]) {
          bdis[tempSnk] = weight;
          bpreLink[tempSnk] = link;
          bQ.push(make_pair(weight, tempSnk));
        }
      }
    }
  }

  if (re) {
    W temp_dis;
    int temp;
    int num = Q.len();
    for (int i = 0; i < num; i++) {
      temp = Q[i].second;

      temp_dis = dis[temp] + bdis[temp];
      if (temp_dis < best_dis) {
        best_dis = temp_dis;
        best_current = temp;
      }
    }
    if (best_dis >= inf) {
#ifdef STATIC_TABLE
      getData(vertex_num, inf, did, false);
#endif
      return false;
    }

    current = best_current;

    while (current != src) {
      path.push_back(preLink[current]);
      graph.findRhs(preLink[current], current, current);
    }
    reverse(path.begin(), path.end());
    current = best_current;
    while (current != snk) {
      path.push_back(bpreLink[current]);
      graph.findRhs(bpreLink[current], current, current);
    }
  }
#ifdef STATIC_TABLE
  getData(vertex_num, inf, did, false);
#endif
  return re;
}

template <typename G, typename WV, typename W>
bool bidijkstra_shortest_path(const G &graph, const WV &NW,
                              const vector<bool> &exclude_nodes,
                              const vector<bool> &exclude_links, const int src,
                              const int snk, vector<int> &path, const W inf) {
  typedef pair<W, int> PII;
  path.clear();
  if (src < 0 || snk < 0) {
    return false;
  }
  if (src == snk) {
    return true;
  }
  int vertex_num = graph.getVertex_num();
  if (src >= vertex_num || snk >= vertex_num) {
    return false;
  }

  vector<int> preLink(vertex_num, -1);
  vector<W> dis(vertex_num, inf);
  vector<char> check(vertex_num, 0);

  LESSOR_T<PII> order;
  size_t j, outDegree, inDegree;
  int link, tempSnk = 0;
  int current;
  W weight;

  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
  dis[src] = 0;
  Q.push(make_pair(0.0, src));

  vector<int> bpreLink(vertex_num, -1);
  vector<W> bdis(vertex_num, inf);
  bdis[snk] = 0;
  Fixed_heap<W, int, LESSOR_T<PII>> bQ(order, vertex_num);
  bQ.push(make_pair(0.0, snk));
  bool re = false;
  W best_dis = inf;
  int best_current = -1;

  while (!Q.empty() && !bQ.empty()) {
    PII p = Q.top();
    current = p.second;
    if (check[current]) {
      re = true;
      break;
    }
    check[current] = 1;
    Q.pop();

    outDegree = graph.getOutDegree(current);
    W current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = graph.getAdj(current, j);
      if (exclude_links[link]) {
        continue;
      }

      graph.findRhs(link, current, tempSnk);
      if (exclude_nodes[tempSnk]) {
        continue;
      }
      weight = current_weight + NW[link];

      if (weight < dis[tempSnk]) {
        dis[tempSnk] = weight;
        preLink[tempSnk] = link;
        Q.push(make_pair(weight, tempSnk));
      }
    }
    p = bQ.top();
    current = p.second;
    if (check[current]) {
      re = true;
      break;
    }
    check[current] = 1;
    bQ.pop();
    inDegree = graph.getInDegree(current);
    current_weight = bdis[current];
    for (j = 0; j < inDegree; j++) {
      link = graph.getReAdj(current, j);
      if (exclude_links[link]) {
        continue;
      }

      graph.findRhs(link, current, tempSnk);
      if (exclude_nodes[tempSnk]) {
        continue;
      }
      weight = current_weight + NW[link];
      if (weight < bdis[tempSnk]) {
        bdis[tempSnk] = weight;
        bpreLink[tempSnk] = link;
        bQ.push(make_pair(weight, tempSnk));
      }
    }
  }
  if (re) {
    W temp_dis;
    int temp;
    int num = Q.len();
    for (int i = 0; i < num; i++) {
      temp = Q[i].second;

      temp_dis = dis[temp] + bdis[temp];
      if (temp_dis < best_dis) {
        best_dis = temp_dis;
        best_current = temp;
      }
    }
    current = best_current;

    while (current != src) {
      path.push_back(preLink[current]);
      graph.findRhs(preLink[current], current, current);
    }
    reverse(path.begin(), path.end());
    current = best_current;
    while (current != snk) {
      path.push_back(bpreLink[current]);
      graph.findRhs(bpreLink[current], current, current);
    }
  }
  return re;
}

namespace inc_ksp {
template <typename W>
struct devote_loc {
  int parent;
  int loc;
  W dis;
  vector<int> exclude_links;

  devote_loc() {
    parent = -1;
    loc = -1;
    dis = numeric_limits<W>::max() / 10;
  }
  devote_loc(const int p, const int l, const W d) {
    parent = p;
    loc = l;
    dis = d;
  }
  bool operator<(const devote_loc &other) const { return dis < other.dis; }
};

template <typename G, typename WV, typename W>
class yen_next_path {
 private:
  const G &graph;
  const WV *weights;
  W inf;
  int src;
  int snk;
  vector<vector<int>> paths;
  vector<bool> exclude_links;
  vector<bool> exclude_nodes;

  vector<bool> temp_exclude_links;
  vector<bool> temp_exclude_nodes;

  unordered_set<int> using_edges;

  priority_queue<devote_loc<W>> Q;
  int last_loc;

 public:
  yen_next_path(const G &graph, const int s, const int t, const WV *ws,
                const vector<bool> &excludeL, W infi_value)
      : graph(graph), weights(ws), inf(infi_value) {
    exclude_nodes.resize(graph.getVertex_num(), false);
    exclude_links = excludeL;
    src = s;
    snk = t;
    last_loc = 0;
  }
  bool next_path(vector<int> &path) {
    path.clear();
    if (src == snk) {
      return true;
    }

    if (paths.empty()) {
      if (!bidijkstra_shortest_path(graph, *weights, exclude_nodes,
                                    exclude_links, src, snk, path, inf)) {
        return false;
      } else {
        paths.push_back(path);
        using_edges.insert(path[0]);
        return true;
      }
    }
    temp_exclude_links = exclude_links;
    temp_exclude_nodes = exclude_nodes;
    int i = 0;
    const vector<int> &last_path = paths.back();
    vector<int> temp_path;
    W pre_dis = 0;
    int temp_src = src;

    for (i = 0; i < last_loc; i++) {
      temp_exclude_nodes[temp_src] = true;
      pre_dis += (*weights)[last_path[i]];
      graph.findRhs(last_path[i], temp_src, temp_src);
    }

    for (; i < (int)last_path.size(); i++) {
      temp_exclude_links[last_path[i]] = true;
      if (i == last_loc) {
        for (unordered_set<int>::iterator it = using_edges.begin();
             it != using_edges.end(); it++) {
          temp_exclude_links[*it] = true;
        }
      }
      if (bidijkstra_shortest_path(graph, *weights, temp_exclude_nodes,
                                   temp_exclude_links, temp_src, snk, temp_path,
                                   inf)) {
        devote_loc<W> temp_loc(paths.size() - 1, i,
                               pre_dis + path_cost(*weights, temp_path, inf));
        temp_loc.exclude_links.push_back(last_path[i]);
        if (i == last_loc) {
          temp_loc.exclude_links.insert(temp_loc.exclude_links.end(),
                                        using_edges.begin(), using_edges.end());
        }
        Q.push(temp_loc);
      }

      temp_exclude_nodes[temp_src] = true;
      pre_dis += (*weights)[last_path[i]];
      graph.findRhs(last_path[i], temp_src, temp_src);
    }
    if (Q.empty()) {
      return false;
    }
    devote_loc<W> temp_loc = Q.top();
    Q.pop();
    last_loc = temp_loc.loc;
    int parent = temp_loc.parent;
    path.clear();
    path.insert(path.end(), paths[parent].begin(),
                paths[parent].begin() + last_loc);
    temp_exclude_links = exclude_links;
    temp_exclude_nodes = exclude_nodes;
    const vector<int> &last_path1 = paths[parent];
    vector<int> temp_path1;

    temp_src = src;
    for (int i = 0; i < last_loc; i++) {
      temp_exclude_nodes[temp_src] = true;
      graph.findRhs(last_path1[i], temp_src, temp_src);
    }

    for (vector<int>::iterator it = temp_loc.exclude_links.begin();
         it != temp_loc.exclude_links.end(); it++) {
      temp_exclude_links[*it] = true;
    }

    bidijkstra_shortest_path(graph, *weights, temp_exclude_nodes,
                             temp_exclude_links, temp_src, snk, temp_path, inf);
    path.insert(path.end(), temp_path.begin(), temp_path.end());
    paths.push_back(path);
    using_edges.clear();
    using_edges.insert(temp_loc.exclude_links.begin(),
                       temp_loc.exclude_links.end());
    using_edges.insert(temp_path.front());

    return true;
  }
};

template <typename G, typename WV, typename W>
class yen_ksp {
 private:
  const G &graph;
  const WV &weights;
  W inf;
  vector<bool> exclude_links;

 public:
  yen_ksp(const G &graph, const WV &ws, const W infi_value)
      : graph(graph), weights(ws), inf(infi_value) {
    exclude_links.resize(graph.getLink_num(), false);
    for (int i = 0; i < graph.getLink_num(); i++) {
      if (weights[i] >= inf) {
        exclude_links[i] = true;
      }
    }
  }

  yen_next_path<G, WV, W> next_path(const int src, const int snk) {
    return yen_next_path<G, WV, W>(graph, src, snk, &weights, exclude_links,
                                   inf);
  }
};
}

namespace maxflow {

template <typename G, typename WV, typename CV, typename W, typename C>

C dijkstra(const G &graph, int src, int snk, const WV &weights,
           const vector<W> &backDis, const CV &caps, vector<C> &flows,
           vector<int> &link_path, const W inf) {
  typedef pair<W, int> PII;

  link_path.clear();

  size_t j, outDegree, inDegree, link, next;
  int current;
  int weight, current_weight;
  int vertex_num = graph.getVertex_num();

  vector<W> dis(vertex_num, inf);
  vector<int> preLink(vertex_num, -1);

  LESSOR_T<PII> order;
  Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
  dis[src] = 0;

  C INF = std::numeric_limits<C>::max() / 10;

  Q.push(backDis[src], src);

  while (!Q.empty()) {
    Q.top(current_weight, current);

    if (current == snk) {
      while (current != src) {
        link_path.push_back(preLink[current]);
        graph.findRhs(preLink[current], current, current);
      }
      reverse(link_path.begin(), link_path.end());
      C re = INF;
      for (vector<int>::const_iterator it = link_path.begin();
           it != link_path.end(); it++) {
        re = min(re, caps[*it]);
      }

      return re;
    }
    Q.pop();

    outDegree = graph.getOutDegree(current);

    for (j = 0; j < outDegree; j++) {
      link = graph.getAdj(current, j);

      if (caps[link] > 0) {
        graph.findRhs(link, current, next);

        weight = current_weight + weights[link];

        if (weight < dis[snk] && weight < dis[next]) {
          preLink[next] = link;
          dis[next] = weight;
          Q.push(weight + backDis[next], next);
        }
      }
    }

    inDegree = graph.getInDegree(current);

    for (j = 0; j < inDegree; j++) {
      link = graph.getReAdj(current, j);
      if (flows[link] > 0) {
        graph.findRhs(link, current, next);
        weight = current_weight - weights[link];
        if (weight < dis[snk] && weight < dis[next]) {
          preLink[next] = link;
          dis[next] = weight;
          Q.push_back(weight + backDis[next], next);
        }
      }
    }
  }

  return 0;
}

template <typename G, typename WV, typename CV, typename W, typename C>
pair<C, W> _minCostMaxFlow(const G &graph, int src, int snk, const WV &weights,
                           const CV &caps, vector<vector<int>> &split_flows,
                           vector<C> &flow_sizes, W infi_value,
                           C max_cap = numeric_limits<C>::max()) {
  int vertex_num = graph.getVertex_num();
  int link_num = graph.getLink_num();
  if (src < 0 || src >= vertex_num) {
    return make_pair<C, W>(0, 0);
  }
  if (snk < 0 || snk >= vertex_num) {
    return make_pair<C, W>(0, 0);
  }
  if (src == snk) {
    return make_pair<C, W>(0, 0);
  }

  vector<W> dist(vertex_num, infi_value);
  vector<C> width(link_num);
  vector<C> flow(link_num);
  vector<W> backDis(vertex_num, infi_value);

  dist[src] = 0;
  width[src] = max_cap;

  dijkstra_shortest_retree(graph, weights, snk, backDis, infi_value);
  CV temp_caps(caps);
  vector<int> link_path;
  C amt = dijkstra(graph, src, snk, weights, backDis, temp_caps, flow,
                   link_path, infi_value);

  while (amt > 0) {
    split_flows.push_back(link_path);
    flow_sizes.push_back(amt);
    for (vector<int>::iterator it = link_path.begin(); it != link_path.end();
         it++) {
      temp_caps[*it] -= amt;
      flow[*it] += amt;
    }

    amt = dijkstra(graph, src, snk, weights, backDis, temp_caps, flow,
                   link_path, infi_value);
  }
}

template <typename G, typename WV, typename CV, typename W, typename C>
pair<C, W> minCostMaxFlow(const G &graph, int src, int snk, const WV &weights,
                          const CV &caps, vector<vector<int>> &split_flows,
                          vector<C> &flow_sizes, W infi_value,
                          C max_cap = numeric_limits<C>::max()) {
  if (graph.isDirect()) {
    return _minCostMaxFlow(graph, src, snk, weights, caps, split_flows,
                           flow_sizes, infi_value, max_cap);
  } else {
    simple_graph temp_graph;
    vector<int> srcs, snks;
    int tempSrc, tempSnk;
    vector<W> tempWights;
    vector<C> tempCaps;
    for (int i = 0; i < graph.getLink_num(); i++) {
      graph.findSrcSnk(i, tempSrc, tempSnk);
      srcs.push_back(tempSrc);
      snks.push_back(tempSnk);

      srcs.push_back(tempSnk);
      snks.push_back(tempSrc);

      tempWights.push_back(weights[i]);
      tempWights.push_back(weights[i]);

      tempCaps.push_back(caps[i]);
      tempCaps.push_back(caps[i]);
    }

    temp_graph.initial(srcs, snks);

    _minCostMaxFlow(temp_graph, src, snk, tempWights, tempCaps, split_flows,
                    flow_sizes, infi_value, max_cap);
  }
}
}
}

#endif
