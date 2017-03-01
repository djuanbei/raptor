/**
 * @file   multicommoditydata.h
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Sun Mar 27 20:40:52 2016
 *
 * @brief  provide data formulti_commodity_flow
 *
 *
 */
#ifndef MULTICOMMODITYDATA_H
#define MULTICOMMODITYDATA_H
#include <limits>
namespace mcmcf {
enum FLOW_TASK_STATUS { UN_ALLOW, SUCCESS_ALLOW, PARTIAL_ALLOW };

enum OPT_STATUS { ALL_OPT, PARTIAL_OPT, ALL_FAILURE };
template <typename G, typename C = double, typename W = double>
class multi_commodity_data {
   public:
    typedef C C_T;
    typedef W W_T;

   private:
    struct flow_task {
        bool isUnsplitted;
        int src, snk;
        C demand;

        flow_task(const int s, const int e, const int d)
            : isUnsplitted(true), src(s), snk(e), demand(d) {}
        void setUnsplitted(bool b) { isUnsplitted = b; }
        bool Unsplitted(void) const { return Unsplitted; }
    };

    struct flow {
        W mean_price;
        vector<vector<size_t>> ver_links;
        vector<C> caps;

        flow() {}
        const vector<size_t>> &getPaths(void) const { return ver_links; }
    };

   private:
    const G &graph;
    const vector<C> &link_capcities;
    const vector<W> &link_prices;

    vector<C> link_current_allows;
    vector<C> link_current_prices;

    W inif;

    Para para;

   public:
    multi_commodity_data(const G &g, const vector<C> &caps,
                         const vector<> &prices)
        : graph(g),
          link_capcities(caps),
          link_prices(prices),
          link_current_allows(caps.size(), 0),
          link_current_prices(prices) {
        inif = numeric_limits<W>::max() / 20;
    }

    double allow_tolerate_rat() const { return 0.3; }
    const G &getTop() const { return graph; }
    C capacity(const int i) const { return link_capcities[i]; }
    W price(const int i) const { return link_prices[i]; }
    C allow(const int i) const { return link_current_allows[i]; }

    const vector<W> &current_price() const { return link_current_allows; }

    W getInif() const { return inif; }
};
}
#endif
