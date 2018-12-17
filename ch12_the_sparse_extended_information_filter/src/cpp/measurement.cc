#include "seif.h"

const static double π = 3.141592653589793238463;

double measure(double);

void measurement(mat& z, std::vector<uint>& c, vec& μ, vec& ξ, mat& Ω, SpMat<ushort>& Λ){
  uint n = (μ.n_elem-3)/2;
  mat qnoise = diagmat(vec({5.0,0.02}));
  Col<uint> cc(c);
  int nnew = max(cc) - n;

  if (nnew > 0) {
    uint osize = ξ.n_elem;
    uint nsize = osize + 2*nnew;
    Ω.resize(nsize,nsize);
    ξ.resize(nsize);
    Λ.resize(n+1+nnew, n+1+nnew);
    μ.resize(nsize);
    mat tmp = inverse_measurement(μ(span(0,2)), z(uvec({0,1}), find(cc > n)));
    μ(span(osize,nsize-1)) = vectorise(tmp);
  }

  cube jcb = jacobian_measurement(μ,c);
  mat zhat = equation_measurement(μ,c);
  mat qnoiseinv = inv(qnoise);

  for (size_t i = 0; i < c.size(); i++) {
    uint j = c[i];
    mat h = jcb.slice(i);
    vec dz = z.col(i)-zhat.col(i);
    dz(1) = measure(dz(1));
    uvec indx = {0,1,2,2*j+1,2*j+2};

    ξ(indx) += h.t()*qnoiseinv*(dz+h*μ(indx));
    Ω(indx,indx) += h.t()*qnoiseinv*h;
  }

  Col<uint> ccc = unique(cc);
  umat edges = zeros<umat>(2,2*ccc.n_elem);
  uint i(0);
  for (size_t j=0; j<ccc.n_elem; ++j, i+=2) {
          edges(1,i) = ccc(j);
          edges(0,i+1) = ccc(j);
          edges(1,i+1) = ccc(j);
  }
  SpMat<ushort> adjacency(edges,ones<Col<ushort>>(2*ccc.n_elem),Λ.n_rows,Λ.n_cols);
  Λ += adjacency + adjacency.t();

}
