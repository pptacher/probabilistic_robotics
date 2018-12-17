#include "seif.h"

void motion(float v, float α, float dt, std::vector<uint> &m0, vec &μ, vec &ξ, mat &Ω, SpMat<ushort> &Λ){
  size_t n = ξ.n_elem;
  mat rnoise = diagmat(vec({8e-4,8e-4,1e-6}));

  mat jcb = jacobian_motion(μ(2),v,α,dt);

  mat ψ = inv(jcb) - eye<mat>(3,3);
  uvec indx = {0,1,2};
  std::vector<uint> m1(m0);
  sort(m1.begin(),m1.end());

  uvec ind(3+2*m1.size());
  ind(uvec({0,1,2})) = uvec({0,1,2});
  uint j(0);
  for (auto iter=m1.begin(); iter != m1.end(); ++iter, j+=2) {
    ind(3+j) = 2**iter+1;
    ind(4+j) = 2**iter+2;
  }
  mat Ω₁ = Ω(ind,uvec(indx));
  mat ψ₁ = Ω₁*ψ;

  mat λ = zeros<mat>(n,n);
  λ(ind,indx) = ψ₁;
  λ(indx,ind) += ψ₁.t();
  λ(indx,indx) += ψ.t()*Ω(indx,indx)*ψ;
  mat Φ = Ω + λ;

  mat Φₓ = Φ(indx,ind);
  mat κ = zeros<mat>(n,n);
  κ(ind,ind) = Φₓ.t()*inv(inv(rnoise)+Φ(indx,indx))*Φₓ;
  Ω += λ - κ;

  vec δ = equation_motion(μ(2),v,α,dt);

  mat lk = λ-κ;
  ξ(ind) += Ω(ind,indx)*δ + lk.rows(ind)*μ;
  μ(span(0,2)) += δ;

  if (!m0.empty()) {
    umat edges = zeros<umat>(2,m1.size()*m1.size());
    uint i(0);
    for (auto iter=m1.begin(); iter != m1.end(); ++iter) {
      for (auto iter1 = m1.begin(); iter1 != m1.end(); ++iter1) {
        edges(0,i) = *iter;
        edges(1,i++) = *iter1;
      }
    }
    if (!edges.empty()) {
      SpMat<ushort> adjacency(edges,ones<Col<ushort>>(edges.n_cols),Λ.n_rows,Λ.n_cols);
      Λ += adjacency;
    }
  }
}
