#include "seif.h"

const static double π = 3.141592653589793238463;
const static float h = 0.74;
const static float l = 2.83;
const static float a = 0.95 + l;
const static float b = 0.50;

mat jacobian_motion(float θ, float v, float α, float dt){
  return {{1, 0, -dt*(v*std::sin(θ)+v/l*std::tan(α)*(a*std::cos(θ)-b*std::sin(θ)))},
          {0, 1,  dt*(v*std::cos(θ)-v/l*std::tan(α)*(a*std::sin(θ)+b*std::cos(θ)))},
          {0, 0,  1}};
}

vec equation_motion(float θ, float v1, float α, float dt){
  float v = v1/(1-std::tan(α)*h/l);
  return {dt*(v*std::cos(θ)-v/l*std::tan(α)*(a*std::sin(θ)+b*std::cos(θ))),
          dt*(v*std::sin(θ)+v/l*std::tan(α)*(a*std::cos(θ)-b*std::sin(θ))),
          dt*v/l*std::tan(α)};
}

cube jacobian_measurement(vec μ, std::vector<uint> indx){
  //pb converting Col<uint> to uvec = Col<long long int>
  uint n = indx.size();
  uvec indx1(n), indx2(n);
  uint i(0);
  for (auto it=indx.begin(); it != indx.end(); ++it, ++i){
    indx1(i) = 2**it+1;
    indx2(i) = 2**it+2;
  }
  vec x = μ(indx1);
  vec y = μ(indx2);

  mat δ = mat(join_horiz(x,y)).each_row()-μ(span(0,1)).t();
  vec q = sum(square(δ),1);

  mat jcb1 = mat({{1,0},
              {0,-1},
              {0, 1},
              {1,0}})*δ.t();
  cube jcb =  zeros<cube>(10,n,1);
  jcb.slice(0) = join_vert(join_vert(join_vert(-jcb1,zeros<rowvec>(n)),-ones<rowvec>(n)),jcb1);
  jcb.slice(0).each_row(uvec({0,2,4,6,8})) %= rowvec(q.t()).for_each( [](mat::elem_type& val){ val = std::sqrt(1/val); });
  jcb.slice(0).each_row(uvec({1,3,7,9})) %= rowvec(q.t()).for_each( [](mat::elem_type& val){ val = 1/val; });

  return reshape(jcb,2,5,n);
}

mat equation_measurement(vec μ, std::vector<uint> indx){
  uint n = indx.size();
  uvec indx1(n), indx2(n);
  uint i(0);
  for (auto it=indx.begin(); it != indx.end(); ++it, ++i){
    indx1(i) = 2**it+1;
    indx2(i) = 2**it+2;
  }
  vec x = μ(indx1);
  vec y = μ(indx2);
  mat δ = mat(join_horiz(x,y)).each_row()-μ(span(0,1)).t();
  vec q = sqrt(sum(square(δ),1));

  mat ϕ = mat(atan2(δ.col(1),δ.col(0))) - μ(2) + π/2;
  return mat(join_vert(q.t(),ϕ.t()));
}

mat inverse_measurement(vec μ, mat z){
  rowvec c = cos(z.row(1)+μ(2));
  rowvec s = sin(z.row(1)+μ(2));

  return mat(join_vert(z.row(0)%s,-z.row(0)%c)).each_col() + μ(span(0,1));
}
