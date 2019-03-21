#include <iostream>
#include <vector>
#include <map>
#include <valarray>
#include <iterator>
#include <algorithm>
#include <limits>
#include <cmath>

typedef unsigned int uint;

struct valarray_comp {
    bool operator()(const std::valarray<int>& lhs, const std::valarray<int>& rhs) const{
      uint s1 = lhs.size();
      uint s2 = rhs.size();
      uint i=0, j=0;
      for ( ;  i < s1 && j < s2; ++i, ++j ) {
        if ( lhs[i] < rhs[j] ) {
          //std::cout << "true" << lhs[0] << ' ' << lhs[1] << ' ' << rhs[0] << ' ' << rhs[1] <<  '\n';
          return true;
        }
        if ( rhs[j] < lhs[i] ) {
          //std::cout << "false" << lhs[0] << ' ' << lhs[1] << ' ' << rhs[0] << ' ' << rhs[1] <<  '\n';
          return false;
      }
    }
    //std::cout << "?" << lhs[0] << ' ' << lhs[1] << ' ' << rhs[0] << ' ' << rhs[1] <<  '\n';
    return ( i==s1) && ( j<s2 );
  }
};

class MarkovDP{

public:

  /*MarkovDP():
    size(10)
    {
      jvalues.assign(10,0);
    }

  MarkovDP(uint n):
    size(n)
    {
      jvalues.assign(n,0);
    }*/

  ~MarkovDP(){}

  void value_iteration(double, uint, std::vector<double>&, uint&) const;
  void add_move(void (*)(const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost));
  void add_vertex(std::valarray<int>,double);

private:
  std::map<std::valarray<int>, double, valarray_comp> vertices;
  std::vector<void (*)(const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&)> moves;

};
