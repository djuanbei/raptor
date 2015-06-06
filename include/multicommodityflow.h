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
#include"composlist.h"
#include<cmath>

using namespace fast_graph;

namespace mcmcf{
  


  enum FLOW_TASK_STATUS{
  
    UN_ALLOW, SUCCESS_ALLOW, 
    PARTIAL_ALLOW
  
  };


  class flow_task{
  private:
    bool isSingle;
    size_t src, snk;
    double demand;
  public:
    flow_task(const size_t  s, const size_t e, const double d ):isSingle( false ), src( s ), snk( e ), demand( d ){
    }
    void setIsSingle(  bool b ){
      isSingle=b;
    }
    bool getIsSingle( void ){
  
      return isSingle;
    }

    friend class multi_commodity_flow;
  
  };

  class flow{

  private:

    double size;
    vector<size_t> ver_links;


  public:
    flow(  ):size( 0 ){
    }

    vector<size_t> getPaths(void  ){
      return ver_links;
    }
    double   getSize( void ) const{
      return size;
    }

    
    friend class multi_commodity_flow;

  };

  enum OPT_STATUS{
    ALL_OPT, PARTIAL_OPT, ALL_FAILURE
  };

  class multi_commodity_flow{

  
  
  private:
    double alpha;
    double desired_cost;
    compressed_sparse_row_graph<int,float,float>& graph;
    double cost;
    size_t link_num;
    vector<flow_task> tasks;
    vector<vector<flow> > flows;
    vector<double> linkCost;
    vector<double> allAllow;


    inline void  differPhi(  ){
      size_t i;
      for (i = 0; i < link_num; i++) {

        allAllow[3*i]=alpha/allAllow[3*i+1] *exp (alpha*allAllow[3*i+2 ]/allAllow[3*i+1 ] );
      }
    }

    double phiValue( void  ) const;
    void initial(  ){
      link_num=graph.getLink_num(  );
      allAllow.resize(2*link_num, 0  );
      linkCost.resize(link_num, 0  );
      //    differPhi.resize(link_num,0 );
      size_t i;
      for (i = 0; i <link_num; i++) {
        linkCost[ i ]=graph.getWeight( i );
        allAllow[ 3*i+1 ]=graph.getCapacity( i );
      }
    }
  

  public:
    multi_commodity_flow(compressed_sparse_row_graph<int,float,float>&g  ):alpha(1),desired_cost( 0 ),  graph( g ), cost( 0 ){
      initial(  );
    
    }

    multi_commodity_flow(const double a, compressed_sparse_row_graph<int,float,float>&g  ):alpha(a),desired_cost( 0 ),  graph( g ), cost( 0 ){
      initial(  );
    }

    multi_commodity_flow(const double a, const double cost, compressed_sparse_row_graph<int,float,float>&g  ):alpha(a),desired_cost( cost ),  graph( g ), cost( 0 ){
      initial(  );

    }

    void setAlpha( const double a ){
  
      alpha=a;
    }
    void setCost(  const double cost ){
      desired_cost=cost;
    }
    double getAlpha( void ) const{
    
      return  alpha;
    }
    double  getCost( void  ) const {
      return desired_cost;
    }


    void addTask(flow_task &task  ){
      tasks.push_back( task );
      vector<flow> dummy;
      flows.push_back( dummy );
    }

    OPT_STATUS optimize(  );

    vector<flow>& getFlowSolutionById( const size_t id ){
      assert( id< flows.size(  ) );
      return flows[ id ];
    }

  };

}

#endif
