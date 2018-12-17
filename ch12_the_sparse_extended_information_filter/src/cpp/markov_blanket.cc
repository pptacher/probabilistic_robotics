#include "seif.h"

std::vector<uint> bfs(uint, SpMat<ushort>&);

std::vector<uint> markov_blanket(uint n,std::vector<uint>& m0, SpMat<ushort>& Λ){
  uint sg = size(Λ,0);
  std::vector<uint> mb1(m0);
  uvec mb = find(Row<ushort>(Λ.row(n)) != 0);

  std::vector<uint> mb2(mb.n_elem);
  uint j=0;
  if (mb(0) != 0) {
    mb2.push_back(0u);
    mb2[j++] = 0;
  }
  for(size_t i=0; i<mb.n_elem; ++i){
    mb2[i+j] = mb(i);
  }

  std::vector<uint> mb12, mb4;
  std::sort(mb1.begin(),mb1.end());
  std::set_intersection(mb1.begin(),mb1.end(),
                        mb2.begin(),mb2.end(),
                        std::back_inserter(mb12));
  if (mb12.empty()) {
    std::vector<uint> mb3;
    mb3 = bfs(n,Λ);
    if (!mb3.empty()) {//should be useless.
      std::sort(mb3.begin(),mb3.end());
      std::set_union(mb3.begin(),mb3.end(),
                      mb2.begin(),mb2.end(),
                      std::back_inserter(mb4));
    }
  }
  else {
    mb4 = mb2;
  }

  std::vector<uint> mb5;
  std::set_union(mb4.begin(),mb4.end(),
                  mb1.begin(),mb1.end(),
                  std::back_inserter(mb5));

  return mb5;
}

std::vector<uint> bfs(uint n, SpMat<ushort>& Λ){
  uint sg = size(Λ,0);
  std::vector<int> vistd(sg,-1);
  vistd[0] = 0;
  std::queue<uint> q;
  q.push(0);
  bool found(false);
  uint currt(0);

  while (!found && !q.empty()) {
    currt = q.front();
    q.pop();
    for (size_t i = 0; i < sg; i++) {
      if (i != currt && Λ.row(currt)(i) && vistd[i]==-1) {
        q.push(i);
        vistd[i] = currt;
      }
    }
    found = vistd[n]==static_cast<int>(currt);
  }

  std::vector<uint> path = {0};
  if (found) {
    path.push_back(currt);
    int prev = vistd[currt];
    while (prev) {
      path.push_back(prev);
      prev = vistd[prev];
    }
  }

  return path;
}
