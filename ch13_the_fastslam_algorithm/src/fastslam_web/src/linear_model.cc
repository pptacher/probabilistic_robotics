#include <armadillo>

using namespace arma;

const static double π = 3.141592653589793238463;
const static float h = 0.74;
const static float l = 2.83;
const static float a = 0.95 + l;
const static float b = 0.50;

double measure(double θ){
  int ϵ = (std::sin(θ) >=0) ? 1 : -1;
  double res =  ϵ*std::acos(std::cos(θ));
  if  (res >= π || res < -π){
    std::cout << ":()" << '\n';
  }
  return res;
}

mat inverse(mat q){
  double d = q(0,0)*q(1,1)-q(0,1)*q(1,0);
  return {{q(1,1)/d,-q(0,1)/d},{-q(1,0)/d,q(0,0)/d}};
}

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

mat jacobian_measurement(vec μ, vec m){
  //pb converting Col<uint> to uvec = Col<long long int>
  vec δ = m-μ(span(0,1));
  vec δ1 = sum(square(δ),0);
  double q = δ1(0);

  mat jcb = {{-δ(0),-δ(1),0,δ(0),δ(1)},
             {δ(1),-δ(0),-q,-δ(1),δ(0)}};
  jcb.row(0) *= (1/sqrt(q));
  jcb.row(1) *= (1/q);

  return jcb;
}

vec equation_measurement(vec μ, vec m){

  vec δ = m-μ(span(0,1));
  vec δ1 = sqrt(sum(square(δ),0));
  double q = δ1(0);
  double ϕ = atan2(δ(1),δ(0)) - μ(2) + π/2;
  return {q,ϕ};
}

mat inverse_measurement(vec μ, mat z){
  rowvec c = cos(z.row(1)+μ(2));
  rowvec s = sin(z.row(1)+μ(2));

  return mat(join_vert(z.row(0)%s,-z.row(0)%c)).each_col() + μ(span(0,1));
}
