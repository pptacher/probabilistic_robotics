#include <cmath>
#include "seif.h"

const static double π = 3.141592653589793238463;
const static double c1 = 1.39;
const static double c2 = 5.99;

SpMat<double> f_xm(uint,uint);
double measure(double);

std::vector<uint> correspondence(mat& z, std::vector<uint>& m0, vec& μ, vec& ξ, mat& Ω, SpMat<ushort>& Λ){
  uint k = size(z,1);

  if (k == 0) {
    return {};
  }

  uint n = (μ.n_elem-3)/2;
  mat qnoise = diagmat(vec({5.0,0.02}));
  Col<uint> c(k);
  uvec newl;
  uvec new1;

  if (n > 0) {
    uint a(0);
    uint m(0);
    std::vector<uint> lst;

    for(auto iter=μ.begin()+3; iter!=μ.end(); iter+=2){
      m++;
      float t1 = *iter-μ(0);
      float t2 = *(iter+1)-μ(1);
      if ( std::abs(t1) < 75 && std::abs(t2) < 75
          && t1*std::cos(μ(2)) + t2*std::sin(μ(2)) > -5 ) {
        ++a;
        lst.push_back(m);
      }
    }

    mat zh = equation_measurement(μ,lst);
    cube jcb = jacobian_measurement(μ,lst);
    mat q = zeros<mat>(2*a,2);
    uint i(0);

    for (auto iter=lst.begin(); iter!=lst.end(); ++iter,++i) {
      std::vector<uint> mkv = markov_blanket(*iter,m0,Λ);
      uvec indx = conv_to<uvec>::from(mkv);

      uvec ind = zeros<uvec>(1+2*indx.n_elem);
      ind(span(0,2)) = uvec({0,1,2});
      for(uint ii=1; ii < indx.n_elem; ++ii){
        ind(span(2*ii+1,2*ii+2)) = uvec({2*indx(ii)+1 ,2*indx(ii)+2});
      }
      uvec ind1 = find(ind(span(3,ind.n_elem-1)) == 2**iter+1,1);
      SpMat<double> f = f_xm(ind1(0),ind.size());

      mat j = jcb.slice(i);
      superlu_opts opts;

      //opts.allow_ugly  = true;
      //opts.symmetric = true;
      //opts.refine = superlu_opts::REF_NONE;
      //mat s(j*f*spsolve(sp_mat(symmatu(Ω(ind,ind))),f.t()*j.t(),"superlu",opts)+qnoise);
      //mat s(j*f*solve(symmatu(Ω(ind,ind)),f.t()*j.t(),solve_opts::fast)+qnoise);

      mat Ωii(trimatu(Ω(ind,ind)));
      mat jf(f.t()*j.t());

      double *Ωii_mem = Ωii.memptr();
      double *jf_mem = jf.memptr();

      MKL_INT info;
      info = LAPACKE_dposv(LAPACK_COL_MAJOR,'U',1+2*indx.n_elem,2,
                    Ωii_mem,1+2*indx.n_elem,
                    jf_mem,ind.size());
      if (info > 0) {
        exit(1);
      }

      //mat sol(jf_mem,1+2*indx.n_elem,2,false,true);
      mat s(j*f*jf+qnoise);

      vec eigval;
      mat eigvec;
      eig_sym(eigval,eigvec,s);
      eigval.for_each( [](mat::elem_type& val){ val = std::sqrt(1/val); });

      q(span(2*i,2*i+1),span::all) = diagmat(eigval)*eigvec.t();
    }

    mat f = zeros<mat>(2*a,k);
    for (size_t j = 0; j < a; j++) {
      f(span(2*j,2*j+1),span::all) =  z.each_col() - zh.col(j);
      f.row(2*j+1).for_each( [](mat::elem_type& val){ val = measure(val); } );
    }

    std::vector<uint> m1(m0);
    std::sort(m1.begin(),m1.end());
    std::vector<uint> ia; //sorted indices of matched elements.
    std::vector<uint> cc;
    set_intersection_indx(lst.begin(), lst.end(),
                                m1.begin(), m1.end(),
                                std::back_inserter(cc),
                                std::back_inserter(ia));

    mat d1 = zeros<mat>(2*a,k);
    mat d2 = zeros<mat>(ia.size(),k);
    int l(0);
    for(auto iter=ia.begin(); iter!=ia.end(); ++iter,++l){
      int j = *iter;
      d1(span(2*j,2*j+1),span::all) =
          q(span(2*j,2*j+1),span::all)*f(span(2*j,2*j+1),span::all);
      d2.row(l) = square(d1.row(2*j)) + square(d1.row(2*j+1));
    }

    urowvec co(k);
    rowvec mins(k);
    for (size_t j = 0; j < k; j++) {
      co(j) = d2.col(j).index_min();
      mins(j) = d2(co(j),j);
    }

    newl = find( mins > c1);
    new1 = find( mins(newl) < c2);

    if (new1.n_elem) {
      urowvec coo = zeros<urowvec>(new1.n_elem);
      for (size_t v = 0; v < new1.n_elem; v++) {
        coo(v) = co(newl(new1(v)));
      }
      urowvec cooo = unique(coo);
      for (int j = cooo.n_elem-1; j >= 0; --j) {
        ia.erase(ia.begin()+ cooo(j));
      }
    }

    Col<uint> tmp(cc);
    c = tmp(co);

    if ((newl.n_elem > 0 && cc.size() < a) || cc.empty()){
      if (cc.empty()) {
        newl = regspace<uvec>(0,k-1);
      }
      std::vector<uint> tmpib(0);
      auto it2 = ia.begin();
      for (size_t j = 0;
          j < a; ++j) {
        if (it2 == ia.end() || *it2 != j ) {
          tmpib.push_back(j);
        }
        else {
          ++it2;
        }
      }

      Col<uint> ib(tmpib);

      l = 0;
      mat d3 = zeros<mat>(tmpib.size(),newl.n_elem);
      for(auto iter=tmpib.begin(); iter!=tmpib.end(); ++iter,++l){
        uint j = *iter;
        d1(uvec({2*j,2*j+1}),newl) =
            q(span(2*j,2*j+1),span::all)*f(uvec({2*j,2*j+1}),newl);
        d3.row(l) = square(d1(uvec({2*j}),newl)) + square(d1(uvec({2*j+1}),newl));
      }

      if (!tmpib.empty()) {
        Col<uint> co1(newl.n_elem);
        rowvec mins1(newl.n_elem);
        for (size_t j = 0; j < newl.n_elem; j++) {
          co1(j) = d3.col(j).index_min();
          c(newl(j)) = lst[ib(co1(j))];
          mins1(j) = d3(co1(j),j);
        }

        new1 = {};
        newl = newl( find(mins1 > c2) );
      }
    }
  } // if (n>0)
  else {
    newl = regspace<uvec>(0,k-1);
    new1 = {};
  }
  if (new1.n_elem > 0) {
    uvec newtmp(newl.n_elem-new1.n_elem);
    uint j(0);
    uint b(0);
    for(size_t i = 0; i < newl.n_elem; ++i){
      if (j < new1.n_elem && i == new1(j)) {
        ++j;
      }
      else {
        newtmp(b++) = newl(i);
      }
    }
    newl = newtmp;
  }
  if (newl.n_elem > 0) {
    uint j(n);
    for(size_t i = 0; i < newl.n_elem; ++i){
      c(newl(i)) = ++j;
    }
  }
  std::vector<uint> corresp(k);
  for(size_t i = 0; i < k; ++i){
    corresp[i] = c(i);
  }
  return corresp;

}

SpMat<double> f_xm(uint p,uint n){
  umat indx = {{0,1,2,3,4},{0,1,2,3+p,4+p}};
  return SpMat<double>(indx,ones<Col<double>>(5),5,n);
}

double measure(double θ){
  double tmp = std::fmod(θ,2*π);
  return tmp - ((int) (tmp/π))*2*π;
}
