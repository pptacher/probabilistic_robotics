#include "markovdp.h"

void move(const std::valarray<int>&, const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&);
void move_rdm(const std::valarray<int>&, const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&);
void warp(const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&);

int main(int argc, char** argv){

  MarkovDP mdp;
  mdp.add_vertex(std::valarray<int>({3,0}), 0.0);
  mdp.add_vertex(std::valarray<int>({2,0}), 0.0);
  mdp.add_vertex(std::valarray<int>({1,0}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,0}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,1}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,2}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,3}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,4}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,5}), 0.0);
  mdp.add_vertex(std::valarray<int>({0,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({-1,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({-2,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({-3,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({-3,7}), 0.0);
  mdp.add_vertex(std::valarray<int>({1,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({2,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({3,6}), 0.0);
  mdp.add_vertex(std::valarray<int>({3,7}), 0.0);
  mdp.add_vertex(std::valarray<int>({-3,8}), -100.0);
  mdp.add_vertex(std::valarray<int>({3,8}), 100.0);

  mdp.add_move([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({0,1}),from,to,cost);});
  mdp.add_move([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({1,0}),from,to,cost);});
  mdp.add_move([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({0,-1}),from,to,cost);});
  mdp.add_move([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({-1,0}),from,to,cost);});
  mdp.add_move(warp);

  std::vector<double> jvalues;
  double epsilon = 1e-3;
  uint max_iter = 1000;
  uint last_iter;
  mdp.value_iteration(epsilon, max_iter, jvalues, last_iter);

  std::cout << "last iter: " << last_iter << '\n';
  for (size_t i = 0; i < jvalues.size(); i++) {
    std::cout << jvalues[i] << '\n';
  }

  return 0;
}

void move(const std::valarray<int>& dxy, const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost) {
  cost = -1.0;
  if (from[1] >= 7) {
    return;
  }
  to[from + dxy] = 1.0;
}

void move_rdm(const std::valarray<int>& dxy, const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost) {
  cost = -1.0;
  if (from[1] >= 7) {
    return;
  }

  to[from + dxy] = 0.9;
  to[from - dxy] = 0.1/3;
  to[from + std::valarray<int>({dxy[1],dxy[0]})] = 0.1/3;
  to[from - std::valarray<int>({dxy[1],dxy[0]})] = 0.1/3;

}


void warp(const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost) {
  cost = 0.0;
  if (from[1] != 7) {
    return;
  }

  to[std::valarray<int>({3,8})] = 0.5;
  to[std::valarray<int>({-3,8})] = 0.5;

}
