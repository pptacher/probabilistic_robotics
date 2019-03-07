#define ARMA_USE_HDF5
//#define ARMA_USE_SUPERLU 1
#define ARMA_64BIT_WORD
#define DARMA_DONT_USE_WRAPPER

#include <iostream>
#include <armadillo>
//#include <string>
//#include <chrono>

#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <queue>
#include <mutex>
#include "mkl_lapacke.h"
//#include "mkl.h"

/////////////////////////////////////////////////////////////////////
//gdb helpers
template<class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}


template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);
template void print_matrix<arma::Col<uint>>(arma::Col<uint>);
template void print_matrix<arma::uvec>(arma::uvec);
template void print_matrix<arma::rowvec>(arma::rowvec);
template void print_matrix<arma::vec>(arma::vec);
/////////////////////////////////////////////////////////////////////


#include <armadillo>
#include <memory>

using namespace arma;

struct NodeRef {
  vec m;
  mat s;
  bool pflag = false;

  NodeRef(){}

  NodeRef(vec from_m,mat from_s,bool from_pflag):
    m(from_m),
    s(from_s),
    pflag(from_pflag){}
};

class BTree{

  class BTreeNode;
  class BTreeNodeI;
  class GaussianNode;

public:
  enum direction : unsigned char {left, right};

  //void initialize(vec,mat);
  bool add_node(uint,vec,mat);
  bool add_node(vec,mat);
  bool rm_node(uint);
  bool set_node(uint,vec,mat);
  bool get_node(uint,vec&,mat&) const;
  uint get_size() const;
  void print(std::ostream&) const;
  void select(vec,std::map<uint, NodeRef>&) const;
  bool dec_pcount(uint);
  void print_tree(std::ostream& os) const;

private:
  uint height = 0;
  uint size = 0;
  static signed char const ptresholdl = -10;
  //static signed char const ptresholdh = 20;
  std::queue<uint> avail;
  //std::map<uint,signed char> pcounts;
  std::vector<signed char> pcounts;
  std::shared_ptr<BTreeNode> root = nullptr;
  std::shared_ptr<BTree::BTreeNode> build_tree(uint,uint,vec&,mat&/*, signed char pc=0*/);
  void preorder(BTreeNode* p,std::ostream&) const;
  void select_rec(BTreeNode*,vec,std::map<uint, NodeRef>&, uint,uint&) const;
  void print_tree_rec(std::ostream& os,BTreeNode* p, uint h) const;

};

class BTree::BTreeNode {
public:
  virtual BTreeNode* next(direction) const= 0;
  virtual std::shared_ptr<BTreeNode> next_s(direction) const = 0;
  virtual void set_next(direction,std::shared_ptr<BTreeNode>) = 0;
  virtual void print(std::ostream&) const = 0;
  virtual bool get_μΣ(vec&,mat&) const = 0;
  //virtual float dec_pcount(uint) = 0;
  //virtual float get_pcount() const = 0;
};

struct BTree::BTreeNodeI : public BTree::BTreeNode{

  BTreeNodeI(std::shared_ptr<BTreeNode> left, std::shared_ptr<BTreeNode> right){
    successors[0] = std::move(left);
    successors[1] = std::move(right);
    //std::cout << "ctor BTreeNodeI" << '\n';
  }

  /*BTreeNodeI( const BTreeNodeI& from){
    successors[0] = from.successors[0];
    successors[1] = from.successors[1];
    std::cout << "cpy ctor BTreeNodeI" << '\n';
  }*/

  BTreeNodeI(){
    std::cout << "dflt ctor BTreeNodeI" << '\n';
  }

  ~BTreeNodeI(){
    //std::cout << "dtor BTreeNodeI" << '\n';
  }

//private:
public:
  std::shared_ptr<BTreeNode> successors[2];

public:
  BTreeNode* next(direction d)const override{
    //std::cout << "BTreeNodeI next" << '\n';
    return successors[d].get();
  }

  virtual std::shared_ptr<BTreeNode> next_s(direction d) const override{
    //std::cout << "BTreeNodeI next_s" << '\n';
    return successors[d];
  }

  void set_next(direction d, std::shared_ptr<BTreeNode> ptr) override{
    successors[d].reset();
    successors[d] = std::move(ptr);
  }

  void print(std::ostream& os) const override{
    os << "left"
            << " get() = " << successors[0].get()
            << " use_count() = " << successors[0].use_count() << '\n' ;
            os << "right"
                    << " get() = " << successors[1].get()
                    << " use_count() = " << successors[1].use_count() << '\n' ;
  }

  virtual bool get_μΣ(vec& m,mat& s) const override {
    return false;
  }

/*  virtual float dec_pcount(uint) override {
    return 0;
  };

  virtual float get_pcount() const override {
    return 0;
  };*/
};

class BTree::GaussianNode  :public BTree::BTreeNode{
public:
  GaussianNode(vec m,mat s/*, float pc = 0*/)
    /*pcount(pc) */{
    this->μ = m;
    this->Σ = s;
    //std::cout << "ctor GaussianNode" << '\n';
    //μ.print("μ:");
    //Σ.print("Σ:");
  }
  ~GaussianNode(){
    //std::cout << "dtor GaussianNode" << '\n';
    //μ.print("μ:");
    //Σ.print("Σ:");
  }
private:
  vec μ = zeros<vec>(2);
  mat Σ = zeros<mat>(2,2);
  //float pcount;
  //mutable std::mutex pcount_mutex;

public:
  BTreeNode* next(direction) const override{
    return nullptr;
  }

  virtual std::shared_ptr<BTreeNode> next_s(direction d) const override{
    return nullptr;
  }

  void print(std::ostream& os) const override{
    μ.print(os,"μ:");
    Σ.print(os,"Σ:");
  }

  virtual bool get_μΣ(vec& m,mat& s) const override {
    m = μ;
    s = Σ;
    return true;
  }

  void set_next(direction d, std::shared_ptr<BTreeNode> ptr) override{};

  /*virtual float dec_pcount(uint nn) override {
    std::lock_guard<std::mutex> guard(pcount_mutex);
    pcount -= 1/((float) nn);
    return pcount;
  }

  virtual float get_pcount() const override {
    std::lock_guard<std::mutex> guard(pcount_mutex);
    return pcount;
  }*/
};
