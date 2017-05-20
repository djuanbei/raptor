/**
 * @file   graph.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Wed Mar  1 21:42:48 2017
 *
 * @brief  smallest graph  class
 *
 *
 */

#ifndef __GRAPH_H
#define __GRAPH_H
#include <algorithm>
#include <cassert>
#include <vector>

namespace raptor {
using std::vector;

class simple_graph {
 private:
  struct startElement {
    int link;
    int src;
    startElement() : link(0), src(0) {}
  };

  struct tempElment {
    int id;
    int src;
    tempElment() : id(0), src(0) {}
    bool operator<(const tempElment &other) const { return src < other.src; }
  };

  bool directed;
  int vertex_num;
  int link_num;

  vector<int> _srcs;
  // out links
  vector<int> outIndex;
  vector<int> link_ends;
  // in links
  vector<int> inIndex;
  vector<startElement> link_starts;

  vector<int> in2outLink_map;  // link_ends  --> orignal link
  vector<int> out2inLink_map;  // orginal link --> link_end

  bool _findSrc(const int link, int &src) const {
    src = vertex_num + 1;
    if (link < 0 || link >= link_num) return false;
    src = _srcs[link];
    return true;
  }
  inline bool _findSnk(const int link, int &snk) const {
    snk = vertex_num + 1;
    if (link >= link_num) return false;
    snk = link_ends[link];
    return true;
  }

  inline bool _findSrcSnk(const int link, int &src, int &snk) const {
    if (!_findSrc(link, src)) return false;
    _findSnk(link, snk);
    return true;
  }

 public:
  void clear() {
    vertex_num = link_num = 0;
    _srcs.clear();
    outIndex.clear();
    link_ends.clear();
    inIndex.clear();
    link_starts.clear();
    in2outLink_map.clear();
    out2inLink_map.clear();
  }

  void initial(vector<int> &srcs, vector<int> &snks, bool direct = false);

  bool isDirect() const { return directed; }

  inline size_t getVertex_num(void) const { return vertex_num; }

  inline int getLink_num(void) const { return link_num; }

  inline int getOutDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return outIndex[vertex + 1] - outIndex[vertex];
  }

  inline int getInDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return inIndex[vertex + 1] - inIndex[vertex];
  }

  int getAdj(int v, int i) const { return in2outLink_map[outIndex[v] + i]; }

  int getReAdj(int v, int i) const {
    return in2outLink_map[link_starts[inIndex[v] + i].link];
  }

  inline bool findSrc(const int link, int &src) const {
    return _findSrc(out2inLink_map[link], src);
  }
  inline bool findSnk(const int link, int &snk) const {
    return _findSnk(out2inLink_map[link], snk);
  }

  bool findRhs(const int link, const int lhs, int &rhs) const {
    int tempSrc, tempSnk;
    if (!findSrcSnk(link, tempSrc, tempSnk)) {
      return false;
    }
    if (lhs == tempSrc) {
      rhs = tempSnk;
      return true;
    }
    if (lhs == tempSnk) {
      rhs = tempSrc;
      return true;
    }

    return false;
  }

  inline bool findSrcSnk(const int link, int &src, int &snk) const {
    return _findSrcSnk(out2inLink_map[link], src, snk);
  }
};
}
#endif
