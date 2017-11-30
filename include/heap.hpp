/**
 * @file   heap.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Mon Apr  6 20:54:06 2015
 *
 * @brief  the min heap used in dijstra shortest path algorithm with the graph
 *only with nonnegtative
 * edge weight
 *
 *
 */

#ifndef _HEAP_H
#define _HEAP_H
#include <utility>  // std::pair, std::make_pair
#include <vector>
using std::vector;
using std::pair;

namespace raptor {

template <class W, class E, class Comp>
class Fixed_heap {
 private:
  vector<pair<W, E>> heap;  // Fixed_heap of Keys
  vector<int> indices;      // Each Key's position (index) in the Fixed_heap
  Comp lt;  // The heap is a minimum-heap with respect to this comparator
  int size;
  vector<E> pass_nodes;

  // Index "traversal" functions
  static inline int left(int i) { return (i << 1) | 1; }
  static inline int right(int i) { return (i + 1) << 1; }
  static inline int parent(int i) { return (i - 1) >> 1; }
  void percolateUp(int i) {
    pair<W, E> x = heap[i];
    int p = parent(i);

    while (i != 0 && lt(x, heap[p])) {
      heap[i] = heap[p];
      indices[heap[i].second] = i;
      i = p;
      p = parent(p);
    }
    heap[i] = x;
    indices[x.second] = i;
  }

  void percolateDown(int i) {
    pair<W, E> x = heap[i];
    while (left(i) < size) {
      int child = right(i) < size && lt(heap[right(i)], heap[left(i)])
                      ? right(i)
                      : left(i);
      if (!lt(heap[child], x)) break;
      heap[i] = heap[child];
      indices[heap[i].second] = i;
      i = child;
    }
    heap[i] = x;
    indices[x.second] = i;
  }

 public:
  Fixed_heap(const Comp &c) : lt(c), size(0) {}
  Fixed_heap(const Comp &c, int cap)
      : heap(cap), indices(cap, -1), lt(c), size(0) {}

  void resize(int cap) {
    heap.resize(cap);
    indices.resize(cap, -1);
    fill(indices.begin(), indices.end(), -1);
    size = 0;
  }

  bool empty() const { return 0 == size; }
  void push(pair<W, E> k) {
    if (-1 == indices[k.second]) {
      pass_nodes.push_back(k.second);
      heap[size] = k;
      indices[k.second] = size;
      percolateUp(size);
      size++;
    } else {
      heap[indices[k.second]].first = k.first;
      percolateUp(indices[k.second]);
    }
  }
  pair<W, E> &top() { return heap[0]; }

  const pair<W, E> &operator[](const int i) const { return heap[i]; }

  void pop() {
    indices[heap[0].second] = -1;
    size--;
    if (size > 0) {
      heap[0] = heap[size];
      indices[heap[0].second] = 0;
      percolateDown(0);
    }
  }

  int len() const { return size; }

  const vector<E> &getPassNodes() const { return pass_nodes; }

  void clear() {
    for (int i = 0; i < size; i++) {
      indices[heap[i].second] = -1;
    }
    size = 0;
    pass_nodes.clear();
  }
};

template <class E, class Comp>
class SFixed_heap {
 private:
  vector<E> heap;       // Fixed_heap of Keys
  vector<int> indices;  // Each Key's position (index) in the Fixed_heap
  const Comp &lt;  // The heap is a minimum-heap with respect to this comparator
  int size;

  // Index "traversal" functions
  static inline int left(int i) { return (i << 1) | 1; }
  static inline int right(int i) { return (i + 1) << 1; }
  static inline int parent(int i) { return (i - 1) >> 1; }
  void percolateUp(int i) {
    E x = heap[i];
    int p = parent(i);

    while (i != 0 && lt[x] < lt[heap[p]]) {
      heap[i] = heap[p];
      indices[heap[i]] = i;
      i = p;
      p = parent(p);
    }
    heap[i] = x;
    indices[x] = i;
  }

  void percolateDown(int i) {
    E x = heap[i];
    while (left(i) < size) {
      int child = right(i) < size && lt[heap[right(i)]] < lt[heap[left(i)]]
                      ? right(i)
                      : left(i);
      if (!(lt[heap[child]] < lt[x])) break;
      heap[i] = heap[child];
      indices[heap[i]] = i;
      i = child;
    }
    heap[i] = x;
    indices[x] = i;
  }

 public:
  SFixed_heap(const Comp &c, int cap)
      : heap(cap), indices(cap, -1), lt(c), size(0) {}

  bool empty() const { return 0 == size; }
  void push(const E k) {
    if (-1 == indices[k]) {
      heap[size] = k;
      indices[k] = size;
      percolateUp(size);
      size++;
    } else {
      percolateUp(indices[k]);
    }
  }
  E &top() { return heap[0]; }

  void pop() {
    indices[heap[0]] = -1;
    size--;
    if (size > 0) {
      heap[0] = heap[size];
      indices[heap[0]] = 0;
      percolateDown(0);
    }
  }

  void clear() {
    size = 0;
    fill(indices.begin(), indices.end(), -1);
  }
};
}
#endif
