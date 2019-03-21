#include "markovdp.h"

void MarkovDP::set_jvalues(vec jval) {
  jvalues = jval;
}
void MarkovDP::set_transitions(std::vector<tst_data> trst) {
  transitions = trst;
}

void MarkovDP::value_iteration(double ϵ, uint max_iter, vec& jval, uvec& policy, uint& last_iter) {
  bool converged = false;
  uint ns = jvalues.n_elem;
  uint i(0);

  while (!converged && (i < max_iter)) {
    mat jvalmat = zeros<mat>(ns,transitions.size());

    for( uint j=0; j < transitions.size(); ++j ) {
      sp_mat const&  tst_proba = transitions[j].tst_proba;
      vec const& tst_cost = transitions[j].tst_cost;

      jvalmat.col(j) = tst_cost + tst_proba*jvalues;
    }

    //to be modified because some actions are unavailable to some states.
    policy = index_max(mat(jvalmat(span::all,span(0,transitions.size()-2))),1);
    jval = arma::max(jvalmat,1);
    /*for (size_t k = 0; k < ns; k++) {
      jval(k) = jvalmat(k,policy(k));
    }*/

    if ( norm(jvalues-jval) < ϵ ) {
      converged = true;
    }
    ++i;
    jvalues.swap(jval);
  }

  jvalues.swap(jval);
  last_iter = i;
}
