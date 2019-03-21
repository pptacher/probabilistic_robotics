#include "markovdp.h"

void MarkovDP::add_move(void (*fptr)(const std::valarray<int>& from, std::map<std::valarray<int>, double, valarray_comp>& to,double& cost)){
  moves.push_back(fptr);
}

void MarkovDP::add_vertex(std::valarray<int> coord,double jvalue){
  vertices.insert(std::pair<std::valarray<int>,double>(coord,jvalue));
}

void MarkovDP::value_iteration(double ϵ, uint max_iter, std::vector<double>& jvalues, uint& last_iter) const{

  uint n = vertices.size();
  uint k = moves.size();
  double* costs = new double[k];
  uint* nb_num = new uint[k*n];

  uint l=0;
  jvalues.resize(n,0);
  std::vector<double> proba;
  std::vector<uint> succ;

  for( auto const& vtx : vertices ) {
    auto const& coord = vtx.first;
    auto const& jvalue = vtx.second;

    int vertex_indx = (uint) std::distance(vertices.begin(),vertices.find(coord));
    jvalues[vertex_indx] = jvalue;

    for (size_t j = 0; j < k; j++) {

      std::map<std::valarray<int>, double, valarray_comp> m;

      moves[j](coord, m, costs[j]);
      double stay = 0.0;
      uint neighb = 0;

      for( auto const& tovtx : m ) {
        auto const&  tocoord = tovtx.first;
        auto const& toprob = tovtx.second;

        auto it = vertices.find(tocoord);
        if (it != vertices.end()) {
          succ.push_back((uint)std::distance(vertices.begin(),it));
          proba.push_back(toprob);
          neighb++;
        }
        else {
          stay += toprob;
        }
      }
      if (stay && (stay < 1)) {
        succ.push_back(vertex_indx);
        proba.push_back(stay);
        neighb++;
      }

      nb_num[l++] = neighb;

    }
  }

  bool converged = false;
  uint iter = 0;
  std::vector<double> jvalues1(n);
  while ( !converged && (converged=true) && ( iter++<max_iter ) ) {

    l = 0;
    for (size_t i = 0; i < n; i++) {

      double max=-60000.0;

      for (size_t j = 0; j < k; j++) {
          if (nb_num[i*k + j] == 0) {
            continue;
          }

          double curr = costs[j];
          for (size_t t = 0; t < nb_num[i*k + j]; t++) {
            curr += proba[l]*jvalues[succ[l]];
            ++l;
          }

          if (curr > max) {
            max = curr;
          }
      }

      jvalues1[i] = (max>-60000.0)?max:jvalues[i];
      if (std::abs(jvalues1[i]-jvalues[i]) > ϵ) {
        converged = false;
      }
    }

    std::swap(jvalues, jvalues1);
  }

  last_iter = iter;

  delete[] nb_num;
  delete[] costs;
}
