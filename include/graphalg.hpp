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
#include <queue>
#include <set>
#include <vector>

#include "heap.h"

using namespace std;

namespace raptor {

template <typename T>
struct LESSOR_T {
    bool operator()(const T &x, const T &y) const { return x.first < y.first; }
};

template <typename WV, typename W>
W path_cost(const WV &NW, const vector<int> &path, W) {
    W re = 0;
    for (vector<int>::const_iterator it = path.begin(); it != path.end();
         it++) {
        re += NW[*it];
    }
    return re;
}

template <typename G, typename WV, typename W>
bool dijkstra_shortest_path(const G &g, const WV &NW, const int src,
                            const int snk, vector<int> &path, const W inf) {
    typedef pair<W, int> PII;
    path.clear();
    if (src < 0 || snk < 0) return false;
    if (src == snk) return true;
    int vertex_num = g.getVertex_num();
    if (src >= vertex_num || snk >= vertex_num) return false;

    vector<int> preLink(vertex_num, -1);
    vector<W> dis(vertex_num, inf);
    LESSOR_T<PII> order;
    size_t j, outDegree;
    int link, next;
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
                g.findRhs(preLink[current], current, current);
            }
            reverse(path.begin(), path.end());
            return true;
        }
        Q.pop();

        outDegree = g.getOutDegree(current);
        W current_weight = p.first;
        for (j = 0; j < outDegree; j++) {
            link = g.getAdj(current, j);

            weight = current_weight + NW[link];

            g.findRhs(link, current, next);

            if (weight < dis[snk] && weight < dis[next]) {
                dis[next] = weight;
                preLink[next] = link;
                Q.push(make_pair(weight, next));
            }
        }
    }
    return false;
}

template <typename G, typename WV, typename W>
void dijkstra_shortest_tree(const G &g, const WV &NW, const int src,
                            vector<int> &preLink, const W inf) {
    typedef pair<W, int> PII;
    int vertex_num = g.getVertex_num();
    preLink.resize(vertex_num);
    fill(preLink.begin(), preLink.end(), -1);

    if (src < 0) return;

    if (src >= vertex_num) return;

    vector<W> dis(vertex_num, inf);
    LESSOR_T<PII> order;
    size_t j, outDegree;
    int link, next;
    int current;
    W weight;

    Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
    dis[src] = 0;
    Q.push(make_pair(0.0, src));
    while (!Q.empty()) {
        PII p = Q.top();
        current = p.second;

        Q.pop();

        outDegree = g.getOutDegree(current);
        W current_weight = p.first;
        for (j = 0; j < outDegree; j++) {
            link = g.getAdj(current, j);

            weight = current_weight + NW[link];

            g.findRhs(link, current, next);

            if (weight < dis[next]) {
                dis[next] = weight;
                preLink[next] = link;
                Q.push(make_pair(weight, next));
            }
        }
    }
}

template <typename G, typename WV, typename W>
void dijkstra_shortest_tree(const G &g, const WV &NW, const int src,
                            vector<int> &preLink, vector<W> &dis, const W inf) {
    typedef pair<W, int> PII;
    int vertex_num = g.getVertex_num();
    preLink.resize(vertex_num);
    dis.resize(vertex_num);
    fill(preLink.begin(), preLink.end(), -1);
    fill(dis.begin(), dis.end(), inf);
    if (src < 0) return;

    if (src >= vertex_num) return;

    LESSOR_T<PII> order;
    size_t j, outDegree;
    int link, next;
    int current;
    W weight;

    Fixed_heap<W, int, LESSOR_T<PII>> Q(order, vertex_num);
    dis[src] = 0;
    Q.push(make_pair(0.0, src));
    while (!Q.empty()) {
        PII p = Q.top();
        current = p.second;

        Q.pop();

        outDegree = g.getOutDegree(current);
        W current_weight = p.first;
        for (j = 0; j < outDegree; j++) {
            link = g.getAdj(current, j);

            weight = current_weight + NW[link];

            g.findRhs(link, current, next);

            if (weight < dis[next]) {
                dis[next] = weight;
                preLink[next] = link;
                Q.push(make_pair(weight, next));
            }
        }
    }
}

template <typename G, typename WV, typename W>
bool bidijkstra_shortest_path(const G &g, const WV &NW, const int src,
                              const int snk, vector<int> &path, const W inf) {
    typedef pair<W, int> PII;
    path.clear();
    if (src < 0 || snk < 0) return false;
    if (src == snk) return true;
    int vertex_num = g.getVertex_num();
    if (src >= vertex_num || snk >= vertex_num) return false;

    vector<int> preLink(vertex_num, -1);
    vector<W> dis(vertex_num, inf);
    vector<bool> check(vertex_num, false);

    LESSOR_T<PII> order;
    size_t j, outDegree, inDegree;
    int link, next = 0;
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
    if (g.isDirect()) {
        while (!Q.empty() && !bQ.empty()) {
            PII p = Q.top();
            current = p.second;
            if (check[current]) {
                re = true;
                break;
            }
            check[current] = true;
            Q.pop();

            outDegree = g.getOutDegree(current);
            W current_weight = p.first;
            for (j = 0; j < outDegree; j++) {
                link = g.getAdj(current, j);

                weight = current_weight + NW[link];

                g.findSnk(link, next);

                if (weight < dis[next]) {
                    dis[next] = weight;
                    preLink[next] = link;
                    Q.push(make_pair(weight, next));
                }
            }
            p = bQ.top();
            current = p.second;
            if (check[current]) {
                re = true;
                break;
            }
            check[current] = true;
            bQ.pop();
            inDegree = g.getInDegree(current);
            current_weight = bdis[current];
            for (j = 0; j < inDegree; j++) {
                link = g.getReAdj(current, j);
                weight = current_weight + NW[link];
                g.findSrc(link, next);
                if (weight < bdis[next]) {
                    bdis[next] = weight;
                    bpreLink[next] = link;
                    bQ.push(make_pair(weight, next));
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
            check[current] = true;
            Q.pop();

            outDegree = g.getOutDegree(current);
            W current_weight = p.first;
            for (j = 0; j < outDegree; j++) {
                link = g.getAdj(current, j);

                weight = current_weight + NW[link];

                g.findRhs(link, current, next);

                if (weight < dis[next]) {
                    dis[next] = weight;
                    preLink[next] = link;
                    Q.push(make_pair(weight, next));
                }
            }
            p = bQ.top();
            current = p.second;
            if (check[current]) {
                re = true;
                break;
            }
            check[current] = true;
            bQ.pop();
            inDegree = g.getInDegree(current);
            current_weight = bdis[current];
            for (j = 0; j < inDegree; j++) {
                link = g.getReAdj(current, j);
                weight = current_weight + NW[link];
                g.findRhs(link, current, next);
                if (weight < bdis[next]) {
                    bdis[next] = weight;
                    bpreLink[next] = link;
                    bQ.push(make_pair(weight, next));
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
            return false;
        }

        current = best_current;

        while (current != src) {
            path.push_back(preLink[current]);
            g.findRhs(preLink[current], current, current);
        }
        reverse(path.begin(), path.end());
        current = best_current;
        while (current != snk) {
            path.push_back(bpreLink[current]);
            g.findRhs(bpreLink[current], current, current);
        }
    }
    return re;
}

template <typename G, typename WV, typename W>
bool bidijkstra_shortest_path(const G &g, const WV &NW,
                              const vector<bool> &exclude_nodes,
                              const vector<bool> &exclude_links, const int src,
                              const int snk, vector<int> &path, const W inf) {
    typedef pair<W, int> PII;
    path.clear();
    if (src < 0 || snk < 0) return false;
    if (src == snk) return true;
    int vertex_num = g.getVertex_num();
    if (src >= vertex_num || snk >= vertex_num) return false;

    vector<int> preLink(vertex_num, -1);
    vector<W> dis(vertex_num, inf);
    vector<bool> check(vertex_num, false);

    LESSOR_T<PII> order;
    size_t j, outDegree, inDegree;
    int link, next = 0;
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
        check[current] = true;
        Q.pop();

        outDegree = g.getOutDegree(current);
        W current_weight = p.first;
        for (j = 0; j < outDegree; j++) {
            link = g.getAdj(current, j);
            if (exclude_links[link]) continue;

            g.findRhs(link, current, next);
            if (exclude_nodes[next]) continue;
            weight = current_weight + NW[link];

            if (weight < dis[next]) {
                dis[next] = weight;
                preLink[next] = link;
                Q.push(make_pair(weight, next));
            }
        }
        p = bQ.top();
        current = p.second;
        if (check[current]) {
            re = true;
            break;
        }
        check[current] = true;
        bQ.pop();
        inDegree = g.getInDegree(current);
        current_weight = bdis[current];
        for (j = 0; j < inDegree; j++) {
            link = g.getReAdj(current, j);
            if (exclude_links[link]) continue;

            g.findRhs(link, current, next);
            if (exclude_nodes[next]) continue;
            weight = current_weight + NW[link];
            if (weight < bdis[next]) {
                bdis[next] = weight;
                bpreLink[next] = link;
                bQ.push(make_pair(weight, next));
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
            g.findRhs(preLink[current], current, current);
        }
        reverse(path.begin(), path.end());
        current = best_current;
        while (current != snk) {
            path.push_back(bpreLink[current]);
            g.findRhs(bpreLink[current], current, current);
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

    set<int> using_edges;

    priority_queue<devote_loc<W>> Q;
    int last_loc;

   public:
    yen_next_path(const G &g, const int s, const int t, const WV *ws,
                  const vector<bool> &excludeL, W infi_value)
        : graph(g), weights(ws), inf(infi_value) {
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

        for (; i < last_path.size(); i++) {
            if (i == last_loc) {
                for (set<int>::iterator it = using_edges.begin();
                     it != using_edges.end(); it++) {
                    temp_exclude_links[*it] = true;
                }
            }
            if (bidijkstra_shortest_path(graph, *weights, temp_exclude_nodes,
                                         temp_exclude_links, temp_src, snk,
                                         temp_path, inf)) {
                devote_loc<W> temp_loc(
                    paths.size() - 1, i,
                    pre_dis + path_cost(*weights, temp_path, inf));

                if (i == last_loc) {
                    temp_loc.exclude_links.insert(temp_loc.exclude_links.end(),
                                                  using_edges.begin(),
                                                  using_edges.end());
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
            graph.findRhs(last_path[i], temp_src, temp_src);
        }

        for (vector<int>::iterator it = temp_loc.exclude_links.begin();
             it != temp_loc.exclude_links.end(); it++) {
            temp_exclude_links[*it] = true;
        }

        bidijkstra_shortest_path(graph, *weights, temp_exclude_nodes,
                                 temp_exclude_links, temp_src, snk, temp_path,
                                 inf);
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
    yen_ksp(const G &g, const WV &ws, const W infi_value)
        : graph(g), weights(ws), inf(infi_value) {
        exclude_links.resize(g.getLink_num(), false);
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
}

#endif
