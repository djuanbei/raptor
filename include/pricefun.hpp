/**
 * @file   pricefun.hpp
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Sun Mar 27 20:38:47 2016
 *
 * @brief  this is the price function for multi_commodity_flow
 *
 *
 */

#ifndef PRICE_FUN_H
#define PRICE_FUN_H
template <typename D>
class price_fun {
 private:
  D &data;

 public:
  /**
   * \phi= sum_i c_i exp( alpha* a_i/c_i )
   *
   */
  double edge_price(int i) const { return edge_price(i, data.allow(i)); }

  double edge_price(const int i, const double allow) const {
    if (data.capacity(i) < 1e-6) {
      return data.getInif();
    }

    if (allow < data.allow_tolerate_rat() * data.capacity(i)) {
      return data.price(i);
    }
    double temp=data.getAlpha(  )) /data.capacity( i );
  return (   data.price( i )*temp *exp(temp*allow );
  }

  bool get_shortest_path(const int src, const int snk,
                         vector<int> &path) const {
    return bidijkstra_shortest_path(data.getTop(), data.current_price(),
                                    data.getInif(), src, snk, path);
  }
};

#endif
