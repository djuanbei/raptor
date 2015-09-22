/**
 * @file   composlist.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Thu Mar 26 20:54:54 2015
 *
 * @brief    compressed sparse row data structure a good data structure for a
 *fixed topology
 *
 *
 */

#ifndef __COMPOSLIST_H
#define __COMPOSLIST_H

#include <cstdio>
#include <vector>
#include <cassert>
#include <algorithm>
#include "heap.h"
#include <queue>
#include <set>
#include <map>
#include <limits>
#include <iostream>
using std::vector;

using namespace std;

namespace fast_graph {

template <class T>
struct LESSOR {
  bool operator()(const T &x, const T &y) const { return x.first < y.first; }
};

template <class W = float, class E = int>
class compressed_sparse_row_graph {
  typedef pair<W, E> PII;
  typedef compressed_sparse_row_graph<W, E> graph_t;
  struct endElement {
    W weight;
    E snk;
    endElement() : weight(0), snk(0) {}
  };

  struct startElement {
    int link;
    E src;
    startElement() : link(0), src(0) {}
  };

  struct precedence {
    E link;
    unsigned char pid;

    precedence() : link(0), pid(0) {}
    precedence &operator=(const precedence &other) {
      link = other.link;
      pid = other.pid;
      return *this;
    }
  };

  struct tempElment {
    int id;
    E src;
  };

  static bool tempElmentCmp(const tempElment &lhs, const tempElment &rhs) {
    return lhs.src < rhs.src;
  }

 private:
  W infi_value;
  size_t ksp_k;

  int vertex_num;
  int link_num;

  vector<E> _srcs;
  // out links
  vector<E> outIndex;
  vector<endElement> link_ends;
  // in links
  vector<E> inIndex;
  vector<startElement> link_starts;

  vector<E> link_map;    // link_ends  --> orignal link
  vector<E> reLink_map;  // orginal link --> link_end

  vector<precedence> shortPaths;  // allpair shortest path

  bool _findSrc(const int link, E &src) const {

    src = vertex_num + 1;
    if (link < 0 || link >= link_num) return false;
    src = _srcs[link];
    return true;

  }
  inline bool _findSnk(const int link, E &snk) const {
    snk = vertex_num + 1;
    if (link >= link_num) return false;
    snk = link_ends[link].snk;
    return true;
  }

  inline bool _findSrcSnk(const int link, E &src, E &snk) const {
    if (!_findSrc(link, src)) return false;
    _findSnk(link, snk);
    return true;
  }

  /**
   *
   *
   * @param src
   * @param snk
   * @param path  the path contain only link
   *
   * @return  true if there have a path connect src and snk
   * false otherwise
   */
  bool _getShortPath(const E src, const E snk,
                     vector<vector<E> > &paths) const {
    assert(src < vertex_num && snk < vertex_num);
    assert(src != snk);
    assert(paths.empty());
    size_t i;
    bool re = false;
    unsigned char pid;
    int link;

    const size_t shift = src * vertex_num;

    for (i = 0; i < ksp_k; i++) {
      if (shortPaths[(shift + snk) * ksp_k + i].link > link_num) return re;

      re = true;
      pid = i;
      E current = snk;
      vector<E> path;
      int num = 0;
      while (num < vertex_num + 2 && current != src) {
        link = (shift + current) * ksp_k + pid;
        path.push_back(shortPaths[link].link);
        _findSrc(shortPaths[link].link, current);
        pid = shortPaths[link].pid;
        num++;
      }
      assert(num < vertex_num + 1);
      std::reverse(path.begin(), path.end());
      paths.push_back(path);
    }

    return true;
  }

 public:
  template <typename T>
  class vertex_map {
   private:
    graph_t &graph;
    vector<T> values;

   public:
    class map_iterator {
     private:
      vertex_map<T> &v_map;
      int index;

     public:
      map_iterator(vertex_map<T> &vv_map, int i = 0)
          : v_map(vv_map), index(i) {}
      map_iterator(map_iterator &other)
          : v_map(other.v_map), index(other.index) {}
      void operator()(map_iterator &other) {
        v_map = other.v_map;
        index = other.index;
      }
      bool operator!=(const map_iterator &other) const {
        return index != other.index;
      }
      void operator++() { index++; }
      void operator--() { index--; }
      void operator++(int) { index++; }

      void operator--(int) { index--; }

      T &operator*() {
        static T error_return;
        if (index < 0 || index >= graph.getVertex_num()) {
          cerr << " there are some thing wrone" << endl;
          return error_return;
        }
        return v_map.values[index];
      }
    };
    typedef map_iterator iterator;

    vertex_map(graph_t &g) : graph(g), values(g.getVertex_num()) {}

    T &operator[](const int index) {
      static T error_return;
      if (index < 0 || index >= graph.getVertex_num()) {
        cerr << " there are some thing wrone" << endl;
        return error_return;
      }
      return values[index];
    }

    const T &operator[](const int index) const {
      static T error_return;
      if (index < 0 || index >= graph.getVertex_num()) {
        cerr << " there are some thing wrone" << endl;
        return error_return;
      }
      return values[index];
    }
    iterator begin() { return map_iterator(*this); }
    iterator end() { return map_iterator(*this, values.size()); }
  };

  template <typename T>
  vertex_map<T> get_vertex_map() {
    return vertex_map<T>(*this);
  }

  compressed_sparse_row_graph(const size_t k = 1)
      : infi_value(numeric_limits<W>::max() / 10e10),
        ksp_k(k),
        vertex_num(0),
        link_num(0) {
    assert(k > 0);
  }
  void clear() {
    vertex_num = link_num = 0;
    _srcs.clear();
    outIndex.clear();
    link_ends.clear();
    inIndex.clear();
    link_starts.clear();
    link_map.clear();
    reLink_map.clear();
    shortPaths.clear();
  }

  void initial(vector<E> &srcs, vector<E> &snks, vector<W> &weights) {
    assert(srcs.size() == snks.size() && srcs.size() == weights.size());
    clear();
    if (0 == srcs.size()) return;
    link_num = srcs.size();
    link_map.resize(link_num);
    reLink_map.resize(link_num);

    vector<tempElment> tContian;

    int i, j, ver_max;
    tempElment dummy_pred;

    ver_max = srcs[0];

    for (i = 0; i < link_num; i++) {
      if (ver_max < srcs[i]) {
        ver_max = srcs[i];
      }
      if (ver_max < snks[i]) {
        ver_max = snks[i];
      }
    }
    vertex_num = ver_max + 1;
    precedence inifPre;
    inifPre.link = link_num + 1;

    shortPaths.resize(ksp_k * vertex_num * vertex_num, inifPre);

    for (size_t i = 0; i < srcs.size(); i++) {
      dummy_pred.id = i;
      dummy_pred.src = srcs[i];
      tContian.push_back(dummy_pred);
    }
    std::sort(tContian.begin(), tContian.end(), tempElmentCmp);

    outIndex.push_back(0);
    i = 0;
    while (i < tContian[0].src) {
      outIndex.push_back(0);
      i++;
    }

    endElement temp;
    temp.weight = weights[tContian[0].id];
    temp.snk = snks[tContian[0].id];

    _srcs.push_back(srcs[tContian[0].id]);
    link_ends.push_back(temp);
    link_map[0] = tContian[0].id;
    reLink_map[tContian[0].id] = 0;

    for (size_t i = 1; i < tContian.size(); i++) {
      if (tContian[i].src != tContian[i - 1].src) {
        for (j = tContian[i - 1].src; j < tContian[i].src; j++) {
          outIndex.push_back(link_ends.size());
        }
      }
      temp.weight = weights[tContian[i].id];
      temp.snk = snks[tContian[i].id];

      _srcs.push_back(srcs[tContian[i].id]);
      link_ends.push_back(temp);
      link_map[i] = tContian[i].id;
      reLink_map[tContian[i].id] = i;
    }

    for (j = tContian[i - 1].src; j < vertex_num; j++) {
      outIndex.push_back(link_ends.size());
    }

    tContian.clear();
    for (size_t i = 0; i < srcs.size(); i++) {
      dummy_pred.id = i;
      dummy_pred.src = snks[i];
      tContian.push_back(dummy_pred);
    }

    std::sort(tContian.begin(), tContian.end(), tempElmentCmp);

    inIndex.push_back(0);
    i = 0;
    while (i < tContian[0].src) {
      inIndex.push_back(0);
      i++;
    }
    startElement dummy;

    dummy.link = reLink_map[tContian[0].id];
    dummy.src = srcs[tContian[0].id];

    link_starts.push_back(dummy);

    for (size_t i = 1; i < tContian.size(); i++) {
      if (tContian[i].src != tContian[i - 1].src) {
        for (j = tContian[i - 1].src; j < tContian[i].src; j++) {
          inIndex.push_back(link_starts.size());
        }
      }

      dummy.link = reLink_map[tContian[i].id];
      dummy.src = srcs[tContian[i].id];
      link_starts.push_back(dummy);
    }

    for (j = tContian[i - 1].src; j < vertex_num; j++) {
      inIndex.push_back(link_starts.size());
    }
  }

  inline size_t getVertex_num(void) const { return vertex_num; }

  void setInfi(const W infi) { infi_value = infi; }

  inline int getLink_num(void) const { return link_num; }

  inline int getOutDegree(E vertex) const {
    assert(vertex < vertex_num);
    return outIndex[vertex + 1] - outIndex[vertex];
  }

  inline int getInDegree(E vertex) const {
    assert(vertex < vertex_num);
    return inIndex[vertex + 1] - inIndex[vertex];
  }
  inline endElement &getLink(E vertex, int k) {
    assert(vertex < vertex_num && k < getOutDegree(vertex));
    return link_ends[outIndex[vertex] + k];
  }

  inline W getWeight(const int link) const {
    return link_ends[reLink_map[link]].weight;
  }

  inline bool findSrc(const int link, E &src) const {
    return _findSrc(reLink_map[link], src);
  }
  inline bool findSnk(const int link, E &snk) const {
    return _findSnk(reLink_map[link], snk);
  }

  inline bool findSrcSnk(const int link, E &src, E &snk) const {
    return _findSrcSnk(reLink_map[link], src, snk);
  }

  inline void setLinkWeight(const int link, const E weight) {
    link_ends[reLink_map[link]].weight = weight;
  }

  void compute_sourceallPair_shortest_path_dijkstra(const E src,
                                                    bool reset = true) {
    const size_t shift = ksp_k * src * vertex_num;

    size_t i, j, k, h, current, outDegree, link, next, next1;
    W weight;

    /**
     * reset the memmory have relation with src
     *
     */

    if (reset) {
      precedence inifPre;
      inifPre.link = link_num + 1;

      fill(shortPaths.begin() + shift,
           shortPaths.begin() + shift + ksp_k * vertex_num, inifPre);
    }

    vector<W> dis(ksp_k * vertex_num, infi_value);
    dis[ksp_k * src] = 0;

    LESSOR<PII> order;

    Fixed_heap<W, E, LESSOR<PII> > Q(order, vertex_num);

    Q.push(make_pair(0.0, src));
    for (i = 0; i < ksp_k; i++) {
      if (i > 0) {
        Q.clear();
        for (int j = 0; j < vertex_num; j++) {
          if (shortPaths[shift + j * ksp_k + i].link < link_num) {
            Q.push(make_pair(dis[j * ksp_k + i], j));
          }
        }
      }

      while (!Q.empty()) {
        PII p = Q.top();
        Q.pop();
        current = p.second;
        outDegree = getOutDegree(current);

        W current_weight = p.first;
        //            shortPaths[ shift+ksp_k*current+i ].weight;

        for (j = 0; j < outDegree; j++) {
          link = outIndex[current] + j;
          const endElement &neighbour = link_ends[link];

          for (k = i; k < ksp_k; k++) {
            next1 = ksp_k * neighbour.snk + k;
            next = shift + next1;
            precedence &dummy_pred = shortPaths[next];
            weight = current_weight + neighbour.weight;
            if (weight < dis[next1]) {
              if (dummy_pred.link < link_num) {
                for (h = ksp_k - 1; h > k; h--) {
                  shortPaths[shift + ksp_k * neighbour.snk + h] =
                      shortPaths[shift + ksp_k * neighbour.snk + h - 1];
                }
              }

              shortPaths[next].link = link;
              shortPaths[next].pid = i;
              dis[next1] = weight;
              if (k == i) Q.push(make_pair(weight, neighbour.snk));
              break;
            }
          }
        }
      }
    }
  }

  /**
   * http://en.wikipedia.org/wiki/Suurballe%27s_algorithm
   *
   * @param src
   * @param snk
   * @param path1  only connnect links
   * @param path2
   */
  void suurballe_shortest(const E src, const E snk, vector<E> &path1,
                          vector<E> &path2) {
    const size_t shift = src * vertex_num;

    /**
     *  fisrt dijkstra shortest path
     *
     */

    compute_sourceallPair_shortest_path_dijkstra(src);

    vector<E> path1link;
    vector<vector<E> > path1links;
    if (!_getShortPath(src, snk, path1links)) {
      return;
    }
    path1link = path1links[0];

    vector<W> newWeight(link_num);
    int i;
    E srcc, snkk;
    /**
     *  new edge weight
     *
     */

    for (i = 0; i < link_num; i++) {
      _findSrcSnk(i, srcc, snkk);
      newWeight[i] = link_ends[i].weight + shortPaths[shift + srcc].weight -
                     shortPaths[shift + snkk].weight;
      assert(newWeight[i] >= 0);
    }

    /**
     *  second dijkstra shortest path
     *
     */

    int link;
    size_t j, current, outDegree;
    W weight;

    vector<W> dis(vertex_num, infi_value);
    vector<E> parent(vertex_num, vertex_num + 10);
    vector<E> shortp(vertex_num, link_num + 10);

    dis[src] = 0;

    LESSOR<PII> order;
    Fixed_heap<W, E, LESSOR<PII> > Q(order, vertex_num);

    Q.push(make_pair(0.0, src));
    while (!Q.empty()) {
      PII p = Q.top();
      current = p.second;
      if (snk == current) {
        break;
      }
      Q.pop();

      outDegree = getOutDegree(current);

      for (j = 0; j < outDegree; j++) {
        link = outIndex[current] + j;
        // reverse all directtion of link in path1
        if (find(path1link.begin(), path1link.end(), link) != path1link.end())
          continue;

        const endElement &neighbour = link_ends[outIndex[current] + j];

        if (src == neighbour.snk) continue;  // cut all link direct to src

        weight = dis[current] + newWeight[link];

        if (weight < dis[neighbour.snk]) {
          shortp[neighbour.snk] = link;
          dis[neighbour.snk] = weight;
          parent[neighbour.snk] = current;

          Q.push(make_pair(weight, neighbour.snk));
        }
      }

      // reverse all directtion of link in path1 and  cut all link direct to src
      for (j = 0; j < path1link.size(); j++) {
        link = path1link[j];
        _findSrcSnk(link, srcc, snkk);
        if (current == snkk && srcc != src) {
          weight = dis[current] + newWeight[link];
          assert(0.0 == newWeight[link]);
          if (weight < dis[srcc]) {
            shortp[srcc] = link;
            dis[srcc] = weight;
            parent[srcc] = current;
            Q.push(make_pair(weight, srcc));
          }
        }
      }

      if (Q.empty()) return;

      vector<E> path2link;
      current = snk;
      while (snk != src) {
        path2link.push_back(shortp[current]);
        current = parent[current];
      }

      /**
       * delete the same links in path1link and path2link
       *
       */

      sort(path2link.begin(), path2link.end());
      sort(path1link.begin(), path1link.end());
      vector<E> sameLink(path1link.size());
      typename vector<E>::iterator it = set_intersection(
          path1link.begin(), path1link.end(), path2link.begin(),
          path2link.end(), sameLink.begin());
      sameLink.resize(it - sameLink.begin());

      for (j = 0; j < sameLink.size(); j++) {
        it = find(path1link.begin(), path1link.end(), sameLink[j]);
        if (it != path1link.end()) {
          path1link.erase(it);
        }
      }

      for (j = 0; j < sameLink.size(); j++) {
        it = find(path2link.begin(), path2link.end(), sameLink[j]);
        if (it != path2link.end()) {
          path2link.erase(it);
        }
      }

      path1link.insert(path1link.end(), path2link.begin(), path2link.end());

      /**
       * obtain the two disjoint paths
       *
       */

      current = src;
      while (current != snk) {
        j = 0;
        for (j = 0; j < path1link.size(); j++) {
          link = path1link[j];
          _findSrc(path1link[j], srcc);
          if (srcc == current) break;
        }

        assert(j < path1link.size());

        path1link.erase(path1link.begin() + j);

        path1.push_back(link);
        _findSnk(link, snkk);
        current = snkk;
      }

      current = src;
      while (current != snk) {
        j = 0;
        for (j = 0; j < path1link.size(); j++) {
          link = path1link[j];
          _findSrc(path1link[j], srcc);
          if (srcc == current) break;
        }

        assert(j < path1link.size());

        path1link.erase(path1link.begin() + j);

        path2.push_back(link);
        _findSnk(link, snkk);
        current = snkk;
      }
    }
  }

  /**
   *
   *  compute two disjoint path1 path2 and the corss srlg do not have
   *intersection
   * @param src
   * @param snk
   * @param srlg map from edge to  srlg
   * @param path1
   * @param path2
   *
   * @return true if find these two paths; false otherwise
   */
  bool twodragonplay(const E src, const E snk, const map<E, vector<int> > &srlg,
                     vector<E> &path1, vector<E> &path2) {
    if (0 == link_num) return true;
    assert(src < vertex_num && snk < vertex_num);
    compute_sourceallPair_shortest_path_dijkstra(src);
    W totalweight = 0.0;
    size_t i;
    for (i = 0; i < link_num; i++) {
      totalweight += getWeight(i);
    }
    double meanweight = (totalweight + 0.0) / link_num;
    map<int, int> srlgoccurnum;

    for (typename map<E, vector<int> >::const_iterator it = srlg.begin();
         it != srlg.end(); it++) {
      for (set<int>::const_iterator iter = it->second.begin();
           it != it->second.end(); it++) {
        if (srlgoccurnum.find(*it) == srlgoccurnum.end()) {
          srlgoccurnum[*it] = 1;
        } else {
          srlgoccurnum[*it]++;
        }
      }
    }

    vector<E> redq, bludq;
    set<int> redsrlg, bluesrlg;
    redq.push_back(src);
    bludq.push_back(src);

    map<int, size_t> redfirstsrlgloc;
    map<int, size_t> bluefirstsrlgloc;

    while (true) {
    }

    return true;
  }



  void compute_allPair_shortest_path() {
    int i;

    precedence inifPre;
    inifPre.link = link_num + 1;

    fill(shortPaths.begin(), shortPaths.end(), inifPre);

#pragma omp parallel for
    for (i = 0; i < vertex_num; i++) {
      compute_sourceallPair_shortest_path_dijkstra(i, false);
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

  bool getShortPath(const E src, const E snk, vector<vector<E> > &paths) const {
    assert(src < vertex_num && snk < vertex_num);
    assert(src != snk);
    assert(paths.empty());

    if (_getShortPath(src, snk, paths)) {
      size_t i, j;
      for (i = 0; i < paths.size(); i++) {
        for (j = 0; j < paths[i].size(); j++) {
          paths[i][j] = link_map[paths[i][j]];
        }
      }

      return true;
    }

    return false;
  }

  bool isValidatePath(const E &src, const E &snk, const vector<E> &path) const {
    E current = src;
    E next = src;
    E next1;
    for (typename vector<E>::const_iterator it = path.begin(); it != path.end();
         it++) {
      if (!findSrcSnk(*it, current, next1)) return false;

      if (current != next) return false;
      next = next1;
    }

    return next == snk;
  }

  static void printPath(const vector<E> &path) {
    if (0 == path.size()) return;

    size_t i = 0;
    string ss;

    while (i < path.size()) {
      std::cout << "by " << path[i] << std::endl;
      i++;
    }
  }
};
}

#endif
