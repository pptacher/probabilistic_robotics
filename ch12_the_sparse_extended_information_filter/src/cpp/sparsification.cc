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
  uint n0(2*m7.size()), n6(2*m8.size());
  uvec indx0(n0), indx6(n6);
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

  MKL_INT info;

  mat Ω₁ = Ω(indx6,
              join_vert(join_vert(uvec({0,1,2}),indx0),indx6));
  mat a1(trimatu(Ω(indx6,indx6)));
  mat Ω₁T(Ω₁.t());

  double *a1_mem = a1.memptr();
  double *b1_mem = Ω₁.memptr();

  info = LAPACKE_dposv(LAPACK_COL_MAJOR,'U',n6,3+n0+n6,
                            a1_mem,n6,
                            b1_mem,n6);
  if (info > 0) {
    exit(1);
  }
  //mat sol1(b1_mem,n6,3+n0+n6,false,true);
  mat l1 = Ω₁T*Ω₁;
  //mat l1 = Ω₁.t()*solve(symmatu(Ω(indx6,indx6)),Ω₁,solve_opts::fast);

  mat Ω₂ = Ω(join_vert(uvec({0,1,2}),indx6),
              join_vert(join_vert(uvec({0,1,2}),indx0),indx6));

  mat Ω₂T(Ω₂.t());
  mat a2(trimatu(Ω(join_vert(uvec({0,1,2}),indx6),join_vert(uvec({0,1,2}),indx6))));

  double *a2_mem = a2.memptr();
  double *b2_mem = Ω₂.memptr();

  info = LAPACKE_dposv(LAPACK_COL_MAJOR,'U',3+n6,3+n0+n6,
                          a2_mem,3+n6,
                          b2_mem,3+n6);
  if (info > 0) {
    exit(1);
  }
  //mat* sol2 = new mat(b2_mem,3+n6,3+n0+n6,false,true);
  mat l2 = Ω₂T*Ω₂;
  //delete sol2;*/
  //mat l2  = Ω₂.t() * solve(symmatu(Ω(join_vert(uvec({0,1,2}),indx6),join_vert(uvec({0,1,2}),indx6))),Ω₂,solve_opts::fast);

  mat Ω₃ = Ω(span(0,2),span::all);
  mat Ω₃T(Ω₃.t());
  mat a3(trimatu(Ω(span(0,2),span(0,2))));

  double *a3_mem = a3.memptr();
  double *b3_mem = Ω₃.memptr();

  info = LAPACKE_dposv(LAPACK_COL_MAJOR,'U',3,Ω.n_rows,
                          a3_mem,3,
                          b3_mem,3);
  if (info > 0) {
    exit(1);
  }
  //mat sol3(b3_mem,3,Ω.n_rows,false,true);
  mat l3 = Ω₃T*Ω₃;

  //mat l3 =Ω₃.t() * solve(symmatu(Ω(span(0,2),span(0,2))),Ω₃,solve_opts::fast);

  mat Ω₄(Ω);
  Ω(join_vert(join_vert(uvec({0,1,2}),indx0),indx6),
    join_vert(join_vert(uvec({0,1,2}),indx0),indx6)) -= l1;
  Ω(join_vert(join_vert(uvec({0,1,2}),indx0),indx6),
    join_vert(join_vert(uvec({0,1,2}),indx0),indx6)) += l2;
  Ω -= l3;

  //correct rounding errors.
  Ω(uvec({0,1,2}),indx6).zeros();
  Ω(indx6,uvec({0,1,2})).zeros();

  ξ += (Ω-Ω₄)*μ;
}
