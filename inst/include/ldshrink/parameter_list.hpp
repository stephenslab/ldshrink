#pragma once

#include <RcppEigen.h>
#include "misc.hpp"

class SingleList{
  using T=int;
  using vectype=std::vector<T>;
  using idtype=Rcpp::StringVector;
  const Rcpp::List input_list;
  const idtype snp_names;
  const  vectype list_a;
  const size_t p;
  const size_t loop_size;
public:
  SingleList(const Rcpp::List  &inp_):input_list(inp_),
				      snp_names(Rcpp::as<idtype>(input_list["snp_id"])),
				      list_a([&](){
					       vectype ret = Rcpp::as<vectype>(input_list["snp_index"]);
					       std::transform(ret.begin(),ret.end(),ret.begin(),
							      [](T t){
								return(t-1);
							      });
					       return(std::move(ret));
					     }()),
				      p(list_a.size()),loop_size((p*p-p)/2+p){}
  vectype get_idx_row() const {
    return (list_a);
  }
  vectype get_idx_col() const {
    return (list_a);
  }
  idtype get_names_col() const {
    return (snp_names);
  }
  idtype get_names_row() const {
    return (snp_names);
  }
  template<int V>
  Rcpp::Vector<V> get_list_elem(const std::string name)const{
    using namespace Rcpp;
    return(as<Vector<V>>(input_list[name]));
  }
  std::pair<T,T> get_idx(const size_t k)const{
    size_t i = floor((2 * p + 1 - sqrt((2 * p + 1) * (2 * p + 1) - 8 * k)) / 2);
    size_t j = k - (2 * p - 1 - i) * i / 2;
    return (std::make_pair(list_a[i], list_a[j]));
  }
  size_t size()const{
    return(loop_size);
  }
};



// class DoubleList{
//   using T=int;
//   using vectype=std::vector<T>;
//   const vectype list_a;
//   const vectype list_b;
//   const size_t p_a;
//   const size_t p_b;
//   const size_t loop_size;
// public:
//   DoubleList(const Rcpp::Vector<cpp2r<T>::data_t> inp_a,const vectype inp_b):
//     list_a(Rcpp::as<vectype>(inp_a)),
//     list_b(Rcpp::as<vectype>(inp_b)),
//     p_a(list_a.size()),
//     p_b(list_b.size()),
//     loop_size(p_a*p_b){
//   }
//   std::pair<vectype,vectype> get_list()const {
//     return(std::make_pair(list_a,list_b));
//   }
//   std::pair<T,T> get_idx(const size_t k)const{
//     size_t i = k %(p_a+1);
//     size_t j = k / p_b;
//     // size_t i = k % p_a
//     // size_t j = k - (2 * p - 1 - i) * i / 2;
//     return (std::make_pair(list_a[i], list_a[j]));
//   }


// };
