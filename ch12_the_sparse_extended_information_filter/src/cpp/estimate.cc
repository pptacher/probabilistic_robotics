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
    }
    μ(span(0,2)) = solve(Ω(span(0,2),span(0,2)),ξ(span(0,2))-Ω(span(0,2),span::all)*μ+Ω(span(0,2),span(0,2))*μ(span(0,2)));
  }
}
