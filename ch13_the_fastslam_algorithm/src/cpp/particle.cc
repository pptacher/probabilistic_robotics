#include <cmath>
#include <iomanip>
#include <random>
#include <thread>
#include "particle.h"
#include "model.h"

const static double π = 3.141592653589793238463;
const static double c2 = 5.99;

void motion_dist(BTree*,mat&,mat&,
                          vec&, mat&,const mat&,
                          const mat&,uint,uint,const uint);

/*void BTree::initialize(vec μ, mat Σ){
  if (root) {
    root.reset();
  }
  uint n = μ.n_elem / 2;
  if (n) {
    size = n;
    height = (uint) ceil(log2(n));
    root = std::move(build_tree(height,(1<<n)-1,(uint) pow(2,height),μ,Σ));
  }
}*/

bool BTree::add_node(vec μ, mat Σ) {
  uint n;
  if (!avail.empty()) {
    n = avail.front();
    avail.pop();
  }
  else {
    n = size;
  }
  return add_node(n,μ,Σ);
}

bool BTree::add_node(uint n,vec μ,mat Σ){
  //std::cout << "adding : " << n << '\n';
  if (root == nullptr) {
    root = std::move(build_tree(0,n,μ,Σ));
    if (size == n) {
      pcounts.push_back(0);
    }
    else
      pcounts[n]=0;
    ++size;
    return true;
  }

  if (n >= (uint) pow(2,height) ) {
    root = std::make_shared<BTree::BTreeNodeI>(
                                  std::move(root),
                                  std::move(build_tree(height,n,μ,Σ)));

    //pcounts[n] = 0;
    /*if (pcounts.size() <= n) {
      pcounts.resize(2*pcounts.size()+1,0);
    }
    pcounts[n]=0;*/
    if (size == n) {
      pcounts.push_back(0);
    }
    else
      pcounts[n]=0;
    ++height;
    ++size;
    return true;
  }

  uint msb = height-1;
  uint m = n >> msb;

  BTreeNode* p = root.get();
  std::shared_ptr<BTreeNode> cpy = std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right));

  BTreeNode* q = cpy.get();

  for (size_t i = 0; i < msb; i++) {

    if (!p->next(static_cast<direction>(m&1))) {
      q->set_next(static_cast<direction>(m&1),std::move(build_tree(msb-i,n,μ,Σ)));
      root = cpy;
      /*if (pcounts.size() <= n) {
        pcounts.resize(2*pcounts.size()+1,0);
      }
      pcounts[n]=0;*/
      if (size == n) {
        pcounts.push_back(0);
      }
      else
        pcounts[n]=0;
      ++size;
      return true;
    }

    p = p->next(static_cast<direction>(m&1));
    q->set_next(static_cast<direction>(m&1),std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right)));
    q = q->next(static_cast<direction>(m&1));
    m = n >> (msb-i-1);
  }

  if (p->next(static_cast<direction>(m&1)) == nullptr) {
    q->set_next(static_cast<direction>(m&1),std::move(build_tree(0,n,μ,Σ)));
    root = cpy;
    /*if (pcounts.size() <= n) {
      pcounts.resize(2*pcounts.size()+1,0);
    }
    pcounts[n]=0;*/
    if (size == n) {
      pcounts.push_back(0);
    }
    else
      pcounts[n]=0;
    ++size;
    return true;
  }

  cpy.reset();
  return false;
}

bool BTree::set_node(uint n,vec μ,mat Σ){

  if (n >= (uint) pow(2,height) || root == nullptr) {
    return false;
  }

  uint msb = height-1;
  uint m = n >> msb;

  BTreeNode* p = root.get();
  std::shared_ptr<BTreeNode> cpy = std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right));

  BTreeNode* q = cpy.get();

  for (size_t i = 0; i < msb; i++) {
    p = p->next(static_cast<direction>(m&1));
    q->set_next(static_cast<direction>(m&1),std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right)));
    q = q->next(static_cast<direction>(m&1));
    m = n >> (msb-i-1);
  }

  if (p->next(static_cast<direction>(m&1)) == nullptr) {
    cpy.reset();
    return false;
  }

  q->set_next(static_cast<direction>(m&1),std::move(build_tree(0,n,μ,Σ)));
  root = cpy;
  /*auto it = pcounts.find(n);
  if (it != pcounts.end()) {
    if (++(it->second) > ptresholdh) {
      pcounts.erase(it);
    }
  }*/
  pcounts[n]++;
  return true;
}

bool BTree::rm_node(uint n) {
  //std::cout << "removing " << n <<  '\n';
  if (n >= (uint) pow(2,height) || root == nullptr) {
    return false;
  }

  uint msb = height-1;
  uint m = n >> msb;

  BTreeNode* p = root.get();
  std::shared_ptr<BTreeNode> cpy = std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right));

  BTreeNode* q = cpy.get();

  for (size_t i = 0; i < msb; i++) {

    p = p->next(static_cast<direction>(m&1));
    q->set_next(static_cast<direction>(m&1),std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right)));
    q = q->next(static_cast<direction>(m&1));
    m = n >> (msb-i-1);

  }

  if (p->next(static_cast<direction>(m&1)) == nullptr) {
    cpy.reset();
    return false;
  }

  q->set_next(static_cast<direction>(m&1),nullptr);
  root = cpy;
  avail.push(n);
  size--;
  return true;
}

bool BTree::get_node(uint n,vec& μ ,mat& Σ) const{

  if (n >= (uint) pow(2,height) || root == nullptr) {
    return false;
  }

  uint msb = height-1;
  uint m = n >> msb;

  BTreeNode* p = root.get();

  for (size_t i = 0; i < msb; i++) {
    p = p->next(static_cast<direction>(m&1));
    m = n >> (msb-i-1);
  }

  if (p->next(static_cast<direction>(m&1)) == nullptr) {
    return false;
  }

  p->next(static_cast<direction>(m&1))->get_μΣ(μ,Σ);
  return true;
}

bool BTree::dec_pcount(uint n) {
  if (n >= (uint) pow(2,height) || root == nullptr) {
    return false;
  }
  //std::cout << pcounts.size() << '\n';
  /*auto it = pcounts.find(n);
  if (it != pcounts.end()) {
    if (--(it->second) < ptresholdl) {
      rm_node(n);
      pcounts.erase(it);
    }
  }*/
    if (--pcounts[n] < ptresholdl) {
    rm_node(n);
    }
    return true;
  /*int msb = height-1;
  uint m = n >> msb;

  BTreeNode* p = root.get();
  std::shared_ptr<BTreeNode> cpy = std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right));

  BTreeNode* q = cpy.get();

  for (size_t i = 0; i < msb; i++) {

    p = p->next(static_cast<direction>(m&1));
    q->set_next(static_cast<direction>(m&1),std::make_shared<BTreeNodeI>(
                                  p->next_s(left),
                                  p->next_s(right)));
    q = q->next(static_cast<direction>(m&1));
    m = n >> (msb-i-1);
  }

  if (p->next(static_cast<direction>(m&1)) == nullptr) {
    cpy.reset();
    return false;
  }

  unsigned char pc =  p->next(static_cast<direction>(m&1))->get_pcount();
  if ( pc-1 < ptresholdl ) {
    q->set_next(static_cast<direction>(m&1),nullptr);
    size--;
    avail.push(n);
  }
  else {
    vec μ(2);
    mat Σ(2,2);
    p->next(static_cast<direction>(m&1))->get_μΣ(μ,Σ);
    q->set_next(static_cast<direction>(m&1),std::move(build_tree(0,n,μ,Σ,pc-1)));
  }
  root = cpy;*/

  /*auto it = pcounts.find(n);
  if (it != pcounts.end()) {
    if (++(it->second) > ptresholdh) {
      pcounts.erase(it);
    }
  }*/
  //pcounts[n]++;

/*  uint msb = height-1;
  uint m = n >> msb;
  uint num_nodes = 1;

  BTreeNode* p = root.get();
  num_nodes *= root.use_count();

  for (size_t i = 0; i < msb; i++) {
    num_nodes *= p->next_s(static_cast<direction>(m&1)).use_count();
    p = p->next(static_cast<direction>(m&1));
    m = n >> (msb-i-1);
  }

  if (p->next(static_cast<direction>(m&1)) == nullptr) {
    return false;
  }

  num_nodes *= p->next_s(static_cast<direction>(m&1)).use_count();
  //decrement persistence count by fractional part of 1 to take into account node sharing between particles.
  double pcount_m = p->next(static_cast<direction>(m&1))->dec_pcount(num_nodes);

  if (pcount_m < (double)ptresholdl) {
    rm_node(n);
  }

  return true;*/
}

std::shared_ptr<BTree::BTreeNode> BTree::build_tree(uint h,uint n,vec& μ,mat& Σ/*, signed char pc*/) {

  if (h == 0) {
    std::shared_ptr<BTree::BTreeNode> res = std::make_shared<BTree::GaussianNode>(μ,Σ);
    return res;
  }

  uint bitmask = 1<<(h-1);
  uint m = bitmask&n;
  std::shared_ptr<BTree::BTreeNode> res;

  if (m) {
    res = std::make_shared<BTree::BTreeNodeI>(
                                  nullptr,
                                  build_tree(h-1,n,μ,Σ));
  }
  else {
    res = std::make_shared<BTree::BTreeNodeI>(
                                  build_tree(h-1,n,μ,Σ),
                                  nullptr);
  }
  return res;
}

void BTree::select(vec μ,std::map<uint, NodeRef>& map) const {
  uint n(0);
  select_rec(root.get(),μ,map,height,n);
}

void BTree::select_rec(BTreeNode* p,vec μ,std::map<uint, NodeRef>& map, uint h,uint& n) const {
  if (!p) {
    n += (1<<h);
    return;
  }
  if (h == 0) {
    vec m(2);
    mat s(2,2);
    p->get_μΣ(m,s);
    vec δ = m - μ(span(0,1));
    if ( std::abs(δ[0]) < 75
        && std::abs(δ[1]) < 75
        && δ[0]*std::cos(μ(2)) + δ[1]*std::sin(μ(2)) > -0.1*norm(δ) ) {
      bool nflag = false;

      if (std::abs(δ[0]) < 30
          && std::abs(δ[1]) < 30
        /*&& δ[0]*std::cos(μ(2)) + δ[1]*std::sin(μ(2)) > 0.01*norm(δ)*/) {
        nflag = true;
      }
      map.try_emplace(n,m,s,nflag);
    }
    ++n;
    return;
  }

  select_rec(p->next(left),μ,map,h-1,n);
  select_rec(p->next(right),μ,map,h-1,n);

}

void BTree::print_tree(std::ostream& os) const {
  BTreeNode* p = root.get();
  print_tree_rec(os,p,height);
}

void BTree::print_tree_rec(std::ostream& os,BTreeNode* p, uint h) const{
  if (!p) {
    return;
  }
  if (h == 0) {
    vec m(2);
    mat s(2,2);
    p->get_μΣ(m,s);
    os << std::setw(15) << m(0)
          << std::setw(15)  << m(1) << '\n';
  }

  print_tree_rec(os,p->next(left),h-1);
  print_tree_rec(os,p->next(right),h-1);
}

void BTree::print(std::ostream& os) const{
  if (!root) {
    return;
  }
  os << "height: " << height
          << " size: " << size << '\n'
          <<"root"
          << " get() = " << root.get()
          << " use_count() = " << root.use_count() << '\n';
  preorder(root.get(),os);
}

void BTree::preorder(BTree::BTreeNode* p,std::ostream& os) const{

  p->print(os);
  os << '\n';
  BTreeNode *leftp = p->next(BTree::left);

  if (leftp) {
    preorder(leftp,os);
  }
  os << '\n';
  BTreeNode*  rightp = p->next(BTree::right);
  if (rightp) {
    preorder(rightp,os);
  }
}

uint BTree::get_size() const {
  return size;
}

std::ostream& operator<< (std::ostream& out, const BTree& t) {
    t.print(out);
    return out;
}

void Particle::motion(double v, double α, double dt, mat z){
  const mat rnoise = diagmat(vec({5e-3,5e-3,2.3e-5}));//2   1e-5 θ   9e-3 x  best 5e-3 2e-5 7e-3&3e-3:NO
  if (z.n_elem == 0) {
    position.each_col( [v,α,dt](vec& a){ a += equation_motion(a[2],v,α,dt);} );
    position += mvnrnd(vec({0,0,0}),rnoise,particleCount);
  }
  else{
    motion_measurement(v,α,dt,z);
    resample();
  }
}

void Particle::motion_measurement(double v, double α, double dt, mat z){
  const mat qnoise = diagmat(vec({0.3,0.015}));//0.8 NO 0.5 OK best 0.3 0.015 O.1:NO 0.005:NO
  const mat rnoise = diagmat(vec({5e-3,5e-3,2.3e-5}));
  mat prediction(position);
  //std::cout << position << '\n';
  prediction.each_col( [v,α,dt](vec& a){ a += equation_motion(a[2],v,α,dt);} );
  //const uint k = z.n_cols;

  const float workLoad = 20.0;
  unsigned int nt = std::min(std::thread::hardware_concurrency(),(uint)std::ceil(particleCount/workLoad));
  std::vector<std::thread> threads;
  for (uint j=0; j<nt; ++j) {
    threads.push_back(std::thread(motion_dist,
                                        b_trees,
                                        std::ref(position),
                                        std::ref(prediction),
                                        std::ref(weight),
                                        std::ref(z),
                                        std::ref(qnoise),
                                        std::ref(rnoise),
                                        j,
                                        nt,
                                        particleCount));
  }

  for (auto& th : threads) {
    th.join();
  }

  threads.clear();

}

void motion_dist(BTree* b_trees,mat& position,mat& prediction,
                          vec& weight, mat& z,const mat& qnoise,
                          const mat& rnoise,uint jj,uint nt,const uint particleCount) {

  uint chunk =  std::ceil(((double)particleCount) / nt);

  std::vector<uint> index;
  for (size_t i = 0; i < particleCount; i+=chunk) {
    index.push_back(i);
  }

  index.push_back(particleCount);
  uint start_indx = index[jj];
  uint end_indx = index[jj+1];
  uint k = z.n_cols;

  /*std::thread::id this_id = std::this_thread::get_id();
    std::cout << "thread " << this_id << " working on"
      <<start_indx << "to" << end_indx <<"\n";*/

  for (uint i = start_indx; i < end_indx; i++){
    BTree* p = &b_trees[i];
    //uint ts = p->get_size();
    //std::cout << "map size: " << ts << '\n';
    vec mini = 1e6*ones<vec>(k);
    vec argmini = zeros<vec>(k);

    std::map<uint, NodeRef> map;
    p->select(prediction.col(i),map);
    uint ms = map.size();
    std::map<uint,uint> indx;
    mat Ω(3*ms,3);
    mat w(ms,k);
    uint ctr = 0;

    for( auto const& [n, ldmrk] : map ) {
      vec μ = ldmrk.m;
      mat Σ = ldmrk.s;

      vec zhat = equation_measurement(prediction.col(i),μ);
      mat h = jacobian_measurement(prediction.col(i),μ);
      mat hx = h(span(0,1),span(0,2));
      mat hm = h(span(0,1),span(3,4));
      mat q = hm*Σ*hm.t()+qnoise;
      mat qinv = inv(q);
      mat Ωx = hx.t()*qinv*hx;
      Ω(span(3*ctr,3*ctr+2),span(0,2)) = Ωx;
      indx[n] = ctr;
      mat δ = z.each_col() - zhat;
      δ.row(1).for_each( [](mat::elem_type& val){ val = measure(val); });
      mat μx = inv(Ωx+inv(rnoise))*hx.t()*qinv*δ;
      vec predx = prediction.col(i);
      mat sample_x(μx);
      sample_x.each_col([Ωx,predx,rnoise](vec& a){a += predx + mvnrnd(a,inv(Ωx+inv(rnoise)),1);});
      mat zhat1(2,k);
      for (size_t l = 0; l < k; l++) {
        zhat1.col(l) = equation_measurement(sample_x.col(l),μ);
      }
      δ = z - zhat1;
      //std::cout << δ << '\n';
      δ.row(1).for_each( [](mat::elem_type& val){ val = measure(val); } );
      for (size_t l = 0; l < k; l++) {
        w(ctr,l) = as_scalar(δ.col(l).t()*qinv*δ.col(l));
      }
      //w.row(j) = δ.each_col([qinv](vec& a){return a.t()*qinv*a;});
      for (size_t l = 0; l < k; l++) {
        if (w(ctr,l) < mini(l)) {
          mini(l) = w(ctr,l);
          argmini(l) = n;
        }
      }
      ++ctr;
    }
    //weight(i) = 1;
    mat Ωx = inv(rnoise);
    vec μx = zeros<vec>(3);
    for (size_t l = 0; l < k; l++) {
      if (mini(l) > c2) {

        //if (z(0,l) > 50 ) {
            //continue;
        //}
        //sample particle position using information gathered so far.
        //vec sample_x =   prediction.col(i) + mvnrnd(μx,inv(Ωx),1);
        vec sample_x = prediction.col(i) + μx;
        vec μ(2);
        mat Σ(2,2);
        μ = inverse_measurement(sample_x,z.col(l));
        mat h = jacobian_measurement(sample_x,μ);
        /*μ = inverse_measurement(prediction.col(i)+μx,z.col(l));
        mat h = jacobian_measurement(prediction.col(i)+μx,μ);*/
        Σ = inv(h(span(0,1),span(3,4))).t()*qnoise*inv(h(span(0,1),span(3,4)));
        mat hm = h(span(0,1),span(3,4));
        mat q = hm*Σ*hm.t()+qnoise;
        weight(i) *= 1/sqrt(det(2*π*q))*exp(-3);
        p->add_node(μ,Σ);
      }
      else  {
        //vec sample_x =   prediction.col(i) + mvnrnd(μx,inv(Ωx),1);
        vec sample_x = prediction.col(i) + μx;
        auto& ldmrk = map[argmini(l)];
        vec μ = ldmrk.m;
        mat Σ = ldmrk.s;
        ldmrk.pflag = false;
        //mat h = jacobian_measurement(prediction.col(i)+μx,μ);
        mat h = jacobian_measurement(sample_x,μ);
        mat hx = h(span(0,1),span(0,2));
        mat hm = h(span(0,1),span(3,4));
        mat q = hm*Σ*hm.t()+qnoise;
        mat qinv = inv(q);
        mat kf = Σ*hm.t()*qinv;
        //vec zhat = equation_measurement(prediction.col(i)+μx,μ);
        vec zhat = equation_measurement(sample_x,μ);
        mat δ = z.col(l) - zhat;
        δ(1) = measure(δ(1));
        μ = μ + kf*δ;
        Σ = (eye<mat>(2,2)-kf*hm)*Σ;
        p->set_node(argmini(l),μ,Σ);
        mat lk = hx*rnoise*hx.t() + hm*Σ*hm.t() + qnoise;
        mat ω = 1/sqrt(det(2*π*lk))*exp(-1/2*δ.t()*inv(lk)*δ);
        weight(i) *= ω(0,0);
        uint j = indx[argmini(l)];
        mat Ωm = Ω(span(3*j,3*j+2),span(0,2));
        Ωx += Ωm;
        μx += inv(Ωx)*hx.t()*qinv*δ;
      }
    }
    //sample particle position using information gathered from measurements.
    position.col(i) = prediction.col(i) + mvnrnd(μx,inv(Ωx),1);

    for( auto const& [n, ldmrk] : map ) {
      if (ldmrk.pflag) {
        p->dec_pcount(n);
      }
    }
  }
}

void Particle::resample(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> v(particleCount);
    for (size_t i = 0; i < particleCount; i++) {
      v[i] = weight(i);
    }
    std::discrete_distribution<> d(v.begin(),v.end());
    std::map<uint, uint> m;
    for(uint n=0; n<particleCount; ++n) {
        ++m[d(gen)];
    }

    uint *cumsum = new uint[particleCount]();
    uint cs = 0;
    uint b=0;
    for( auto const& [k,val] : m ) {
      while (b < k) {
        cumsum[b++] = cs;
      }
      cs += val;
      cumsum[k] = cs;
      ++b;
    }
    while (b < particleCount) {
      cumsum[b++] = cs;
    }
    const float workLoad = 50.0;
    unsigned int nt = std::min(std::thread::hardware_concurrency(),(uint)std::ceil(particleCount/workLoad));
    uint chunk =  std::ceil(((double)particleCount) / nt);

    std::vector<uint> index;
    for (size_t i = 0; i < particleCount; i+=chunk) {
      index.push_back(i);
    }
    index.push_back(particleCount);

    BTree* b_trees_tmp = new BTree[particleCount];
    mat position_tmp(3,particleCount);

    std::vector<std::thread> threads;
    for (uint j=0; j<nt; ++j) {
      threads.push_back(std::thread(&Particle::resample_dist,
                                          this,
                                          b_trees_tmp,
                                          std::ref(position_tmp),
                                          std::ref(m),
                                          cumsum,
                                          index[j],
                                          index[j+1]));
    }

    for (auto& th : threads) {
      th.join();
    }

    threads.clear();

    std::swap(b_trees,b_trees_tmp);
    position = position_tmp;
    weight.ones();
    delete[] b_trees_tmp;
    delete[] cumsum;

}

void Particle::resample_dist(BTree* bt,
                              mat& pt,
                              std::map<uint,uint>& m,
                              uint* cumsum,
                              uint start_indx,
                              uint end_indx) {
  uint i;
  for (i = 0; cumsum[i] <= start_indx ; i++) {}
  uint start_k = i;
  for (;i < particleCount && cumsum[i] <= end_indx ; i++) {}
  uint end_k = i;
  //uint j=start_indx;
  /*std::thread::id this_id = std::this_thread::get_id();
    std::cout << "thread " << this_id << " working on"
      <<start_k << "to" << end_k <<"\n";*/
   for( uint t=start_k; t<end_k; ++t) {
      for (uint j=(t==0)?0:cumsum[t-1];
            j < cumsum[t];
            ++j) {

        bt[j] = b_trees[t];
        pt.col(j) = position.col(t);
      }
    }

}

void Particle::print(std::ostream& os) const{
  for (size_t i = 0; i < particleCount; i++) {
    os << std::setw(15) << position(0,i)
          << std::setw(15)  << position(1,i)
          << std::setw(15)  << position(2,i);
  }
}

std::ostream& operator<< (std::ostream& out, const Particle& t) {
    t.print(out);
    return out;
}
