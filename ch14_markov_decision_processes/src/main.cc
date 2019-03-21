#include "markovdp.h"

void move(const std::valarray<int>&, const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&);
void move_rdm(const std::valarray<int>&, const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&);
void warp(const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&);
std::vector<tst_data> build_transitions(const mat&,
                          std::vector<void (*)(const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&)>&);


int main(int argc, char** argv){

  mat occupancy;
  occupancy.load(arma::hdf5_name("occupancy_map.h5", "occupancy"));
  uint n = 0;
  for (size_t i = 0; i < occupancy.n_rows; i++) {
    for (size_t j = 0; j < occupancy.n_cols; j++) {
      if (occupancy(i,j) == 0) {
        occupancy(i,j) =n++;
        continue;
      }
    }
  }

  vec jvalues = zeros<vec>(54);
  //terminal states.
  jvalues(0) = +100;
  jvalues(1) = -100;
  jvalues(36) = -100;
  jvalues(37) = +100;

  std::vector<void (*)(const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&)> moves;

  moves.push_back([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({0,1,0}),from,to,cost);});
  moves.push_back([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({1,0,0}),from,to,cost);});
  moves.push_back([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({0,-1,0}),from,to,cost);});
  moves.push_back([](const auto& from, auto& to, double& cost){
                  return move_rdm(std::valarray<int>({-1,0,0}),from,to,cost);});
  moves.push_back(warp);

  std::vector<tst_data> transitions = build_transitions(occupancy, moves);

  MarkovDP mdp(jvalues, transitions);

  uvec policy(54);
  uint last_iter(0);
  double epsilon = 1e-2;
  uint max_iter = 1000;
  mdp.value_iteration(epsilon, max_iter, jvalues, policy, last_iter);

  std::cout << "converged in " << last_iter << " iterations." << '\n';
  jvalues.print("Jvalues: ");
  policy.print("optimal policy: ");

  return 0;
}

void move(const std::valarray<int>& dxy, const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost) {
  cost = -1.0;
  if (from[0] == 0 || ( from[0] == 7 && from[1] == 6 && from[2] == 0 )) {
    return;
  }
  to[from + dxy] = 1.0;
}

void move_rdm(const std::valarray<int>& dxy, const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost) {
  cost = -1.0;
  if (from[0] == 0 || ( from[0] == 7 && from[1] == 6 && from[2] == 0 )) {
    return;
  }

  to[from + dxy] = 0.9;
  to[from - dxy] = 0.1/3;
  to[from + std::valarray<int>({dxy[1],dxy[0],dxy[2]})] = 0.1/3;
  to[from - std::valarray<int>({dxy[1],dxy[0],dxy[2]})] = 0.1/3;

}

void warp(const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost) {
  cost = 0.0;

  if (from[0] == 0 && from[2] == 0) {
    to[std::valarray<int>({from[0], from[1], 1})] = 0.5;
    to[std::valarray<int>({from[0], from[1], -1})] = 0.5;
    return;
  }

  if ( from[0] == 7 && from[1] == 6 && from[2] == 0) {
    to[std::valarray<int>({7, 6, 1})] = 0.5;
    to[std::valarray<int>({7, 6, -1})] = 0.5;
    return;
  }

}

inline uint to_linear(std::valarray<int> coord, const mat& occupancy) {
  return (coord[2]+1)*18 + (uint)occupancy(coord[0],coord[1]);
}

std::vector<tst_data> build_transitions(const mat& occupancy,
                          std::vector<void (*)(const std::valarray<int>&, std::map<std::valarray<int>, double, valarray_comp>&,double&)>& moves) {

  std::vector<tst_data> transitions;
  for (size_t i = 0; i < moves.size(); i++) {
    transitions.push_back({sp_mat(54,54),vec(54)});
  }
  int nr = occupancy.n_rows;
  int nc = occupancy.n_cols;

  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      if (occupancy(i,j)==-1) {
        continue;
      }
      for (int z = -1; z < 2; z++) {
        std::valarray<int> from({i, j, z});
        uint lin_indx = to_linear(from, occupancy);
        for ( uint m=0; m < moves.size(); ++m) {

          std::map<std::valarray<int>, double, valarray_comp> to;
          double cost;
          moves[m](from,to,cost);

          double stay(0);
          for( auto const& tovtx : to ) {
            auto const&  tocoord = tovtx.first;
            auto const& toprob = tovtx.second;


            if (tocoord[0] < 0 || tocoord[0] >= nr ||
                tocoord[1] < 0 || tocoord[1] >= nc ||
                occupancy(tocoord[0], tocoord[1]) == -1) {
                stay += toprob;
            }
            else {
              uint lin_indx1 = to_linear(tocoord, occupancy);
              transitions[m].tst_proba(lin_indx,lin_indx1) = toprob;
            }
          }

          transitions[m].tst_proba(lin_indx,lin_indx) = (to.size())?stay:1;
          transitions[m].tst_cost(lin_indx) = (stay<1)?cost:0;

        }
      }
    }
  }

  return transitions;
}
