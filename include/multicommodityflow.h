/**
 * @file   multicommodityflow.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Thu Mar 26 21:02:41 2015
 *
 * @brief  @see http://en.wikipedia.org/wiki/Multi-commodity_flow_problem
 *
 *
 */

#ifndef MULTICOMMODITYFLOW_H
#define MULTICOMMODITYFLOW_H
#include <cmath>
#include "graph.h"

// using namespace fast_graph;

namespace mcmcf {

class multi_commodity_flow {
    enum INDEX {
        CURRENT_PRICE = 0,
        ORIGINAL_PRICE = 1,
        CURRENT_ALLOW = 2,
        ORIGINAL_CAPACITY = 3
    } private :

        double alpha;
    double desired_cost;

    double cost;
    size_t link_num;
    vector<flow_task> tasks;
    vector<vector<flow>> flows;
    vector<double> linkCost;
    /**
     * 0 current link price
     * 1 orignal link price
     * 2 current link allow
     * 3  orignal link capacity
     */

    vector<double> allAllow;

    /**
     * \phi= sum_i c_i exp( alpha* a_i/c_i )
     *
     */
    inline void differPhi() {
        size_t i;
        for (i = 0; i < link_num; i++) {
            size_t index = 4 * i;
            allAllow[index] = alpha / allAllow[index + ORIGINAL_CAPACITY] *
                              allAllow[index + ORIGINAL_PRICE] exp(
                                  alpha * allAllow[index + CURRENT_ALLOW] /
                                  allAllow[index + ORIGINAL_CAPACITY]);
        }
    }

    double phiValue(void) const;
    void initial() {
        link_num = graph.getLink_num();
        allAllow.resize(4 * link_num, 0);
        linkCost.resize(link_num, 0);
        //    differPhi.resize(link_num,0 );
        size_t i;
        for (i = 0; i < link_num; i++) {
            linkCost[i] = graph.getWeight(i);
        }
    }

   public:
    multi_commodity_flow(const G &g)
        : alpha(1), desired_cost(0), graph(g), cost(0) {
        initial();
    }

    multi_commodity_flow(const double a, const G &g)
        : alpha(a), desired_cost(0), graph(g), cost(0) {
        initial();
    }

    multi_commodity_flow(const double a, const double cost, const G &g)
        : alpha(a), desired_cost(cost), graph(g), cost(0) {
        initial();
    }

    void setAlpha(const double a) { alpha = a; }
    void setCost(const double cost) { desired_cost = cost; }
    double getAlpha(void) const { return alpha; }
    double getCost(void) const { return desired_cost; }

    void addTask(flow_task &task) {
        tasks.push_back(task);
        vector<flow> dummy;
        flows.push_back(dummy);
    }

    OPT_STATUS optimize();

    vector<flow> &getFlowSolutionById(const size_t id) {
        assert(id < flows.size());
        return flows[id];
    }
};
}

#endif
