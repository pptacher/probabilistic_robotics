#define ARMA_USE_HDF5
//#define ARMA_USE_SUPERLU 1
#define ARMA_64BIT_WORD
#define DARMA_DONT_USE_WRAPPER

#include <iostream>
#include <armadillo>
//#include <string>

#include <iostream>
#include <vector>
#include <map>
#include <valarray>
#include <iterator>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iomanip>
#include "mkl_lapacke.h"
//#include "mkl.h"

/////////////////////////////////////////////////////////////////////
//gdb helpers
template<class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}


template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);
template void print_matrix<arma::Col<uint>>(arma::Col<uint>);
template void print_matrix<arma::uvec>(arma::uvec);
template void print_matrix<arma::rowvec>(arma::rowvec);
template void print_matrix<arma::vec>(arma::vec);
/////////////////////////////////////////////////////////////////////

using namespace arma;

typedef unsigned int uint;

struct valarray_comp {
    bool operator()(const std::valarray<int>& lhs, const std::valarray<int>& rhs) const{
      uint s1 = lhs.size();
      uint s2 = rhs.size();
      uint i=0, j=0;
      for ( ;  i < s1 && j < s2; ++i, ++j ) {
        if ( lhs[i] < rhs[j] ) {
          return true;
        }
        if ( rhs[j] < lhs[i] ) {
          return false;
      }
    }
    return ( i==s1) && ( j<s2 );
  }
};

struct  tst_data {
  sp_mat tst_proba;
  vec tst_cost;
};

class MarkovDP{

public:

  MarkovDP(){}

  MarkovDP(vec jval, std::vector<tst_data> tst):
    jvalues(jval),
    transitions(tst){}

  ~MarkovDP(){}

  void value_iteration(double, uint, vec&, uvec&, uint&);
  void set_jvalues(vec);
  void set_transitions(std::vector<tst_data>);

private:
  vec jvalues;
  std::vector<tst_data> transitions;

};
