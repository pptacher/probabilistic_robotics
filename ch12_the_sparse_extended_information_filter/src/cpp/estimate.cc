#include "seif.h"

void estimate(std::vector<uint> &m0, vec &μ, vec &ξ, mat &Ω){
  uint n = (μ.n_elem-3)/2;

  if (n < 100) {
    μ = solve(Ω,ξ);
  }
  else {
    for (auto iter = m0.begin(); iter != m0.end(); ++iter) {

      auto indx = span(2**iter+1,2**iter+2);
      μ(indx) = solve(Ω(indx,indx),ξ(indx)-Ω(indx,span::all)*μ+Ω(indx,indx)*μ(indx));

      /*mat a(trimatu(Ω(indx,indx)));
      mat b(ξ(indx)-Ω(indx,span::all)*μ+Ω(indx,indx)*μ(indx));

      double *a_mem = a.memptr();
      double *b_mem = b.memptr();

      MKL_INT info;
      info = LAPACKE_dposv(LAPACK_COL_MAJOR,'U',2,1,
                    a_mem,2,
                    b_mem,2);
      if (info > 0) {
        exit(1);
      }
      μ(indx) = mat(b_mem,2,1,true);*/
    }
    μ(span(0,2)) = solve(Ω(span(0,2),span(0,2)),ξ(span(0,2))-Ω(span(0,2),span::all)*μ+Ω(span(0,2),span(0,2))*μ(span(0,2)));
  }
}
