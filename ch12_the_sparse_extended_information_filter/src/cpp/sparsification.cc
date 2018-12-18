#include "seif.h"

void sparsification(std::vector<uint>& m0, std::vector<uint>& m6, vec& μ, vec& ξ, mat& Ω, SpMat<ushort>& Λ){

  std::vector<uint> m7(m0);
  std::vector<uint> m8(m6);
  std::vector<uint> m9;

  std::sort(m7.begin(), m7.end());
  std::sort(m8.begin(), m8.end());
  std::set_union(m7.begin(), m7.end(),
                  m8.begin(), m8.end(),
                  std::back_inserter(m9));

  for (auto it1 = m6.begin(); it1 != m6.end(); ++it1 ){
    Λ(0,*it1) = 0;
    Λ(*it1,0) = 0;
    for (auto it2 = m9.begin(); it2 != m9.end(); ++it2){
      if (*it1 != *it2) {
        Λ(*it1,*it2) = 1;
        Λ(*it2,*it1) = 1;
      }
    }
  }
  uvec indx0(2*m7.size()), indx6(2*m8.size());
  uint i(0);
  for (auto it=m7.begin(); it != m7.end(); ++it,i+=2){
    indx0(i) = 2**it+1;
    indx0(i+1) = 2**it+2;
  }
  i=0;
  for (auto it=m8.begin(); it != m8.end(); ++it, i+=2){
    indx6(i) = 2**it+1;
    indx6(i+1) = 2**it+2;
  }

  mat Ω₁ = Ω(join_vert(join_vert(uvec({0,1,2}),indx0),indx6),
              indx6);
  mat l1 = Ω₁*solve(symmatu(Ω(indx6,indx6)),Ω₁.t(),solve_opts::fast);

  mat Ω₂ = Ω(join_vert(join_vert(uvec({0,1,2}),indx0),indx6),
              join_vert(uvec({0,1,2}),indx6));
  mat l2  = Ω₂ * solve(symmatu(Ω(join_vert(uvec({0,1,2}),indx6),join_vert(uvec({0,1,2}),indx6))),Ω₂.t(),solve_opts::fast);

  mat Ω₃ = Ω(span::all,span(0,2));
  mat l3 =Ω₃ * solve(symmatu(Ω(span(0,2),span(0,2))),Ω₃.t(),solve_opts::fast);

  mat Ω₄(Ω);
  Ω(join_vert(join_vert(uvec({0,1,2}),indx0),indx6),
    join_vert(join_vert(uvec({0,1,2}),indx0),indx6)) -= l1;
  Ω(join_vert(join_vert(uvec({0,1,2}),indx0),indx6),
    join_vert(join_vert(uvec({0,1,2}),indx0),indx6)) += l2;
  Ω -= l3;

  Ω(uvec({0,1,2}),indx6).zeros();
  Ω(indx6,uvec({0,1,2})).zeros();

  ξ += (Ω-Ω₄)*μ;
}
