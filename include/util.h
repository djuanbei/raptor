#include <algorithm>
#include<vector>
namespace raptor {
using namespace std;

/** 
 * 
 * 
 * @param first Forward iterators to the initial and final positions of a sorted (or properly partitioned) sequence. The range used is [first,last), which contains all the elements between first and last, including the element pointed by first but not the element pointed by last.
 * @param last  elements in [first, last) must be sorted 
 * @param value 
 * 
 * @return non-negative integer if value equal to one element in [first , last)
 */
template<typename T>
int binfind(typename vector<T>::const_iterator first, typename vector<T>::const_iterator last, T  value){
  typename vector<T>::const_iterator it=lower_bound(first, last, value);
  if(it==last){
    return -1;
  }
  if(*it!= value){
    return -1;
  }
  return it-first;
  
}

template<typename T>
int binfind(typename vector<T>::iterator first, typename vector<T>::iterator last, T  value){
  typename vector<T>::iterator it=lower_bound(first, last, value);
  if(it==last){
    return -1;
  }
  if(*it!= value){
    return -1;
  }
  return it-first;
  
}

}
