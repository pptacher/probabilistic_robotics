#define ARMA_USE_HDF5
//#define ARMA_USE_SUPERLU 1
#define ARMA_64BIT_WORD
#define DARMA_DONT_USE_WRAPPER

#include <iostream>
#include <armadillo>
#include <string>
#include <chrono>
#include <boost/progress.hpp>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <queue>
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

void seif(void);
void motion(float, float, float, std::vector<uint>&, vec&, vec&, mat&, SpMat<ushort>&);
void estimate(std::vector<uint>&, vec&, vec&, mat&);
std::vector<uint> correspondence(mat&, std::vector<uint>&, vec&, vec&, mat&, SpMat<ushort>&);
void measurement(mat&, std::vector<uint>&, vec&, vec&, mat&, SpMat<ushort>&);
void sparsification(std::vector<uint>&, std::vector<uint>& ,vec&, vec&, mat&, SpMat<ushort>&);

std::vector<uint> markov_blanket(uint, std::vector<uint>&, SpMat<ushort>&);

mat jacobian_motion(float, float, float, float);
vec equation_motion(float, float, float, float);
cube jacobian_measurement(vec, std::vector<uint>);
mat equation_measurement(vec, std::vector<uint>);
mat inverse_measurement(vec, mat);

template<class InputIt1, class InputIt2, class OutputIt1, class OutputIt2>
OutputIt1 set_intersection_indx(InputIt1 first1, InputIt1 last1,
                          InputIt2 first2, InputIt2 last2,
                          OutputIt1 d_first, OutputIt2 indx){
    auto tmp = first1;
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            ++first1;
        } else  {
            if (!(*first2 < *first1)) {
                *indx++ = std::distance(tmp, first1);
                *d_first++ = *first1++;
            }
            ++first2;
        }
    }
    return d_first;
}
