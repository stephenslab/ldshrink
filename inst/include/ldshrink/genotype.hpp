#pragma once
#include "misc.hpp"
#include <RcppEigen.h>
#include <memory>
#include <Rinternals.h>

class Genotype{
  const int RTYPE= INTSXP;
  using T=int;
  using data_v = std::pair<const Eigen::VectorXd,double>;
  using	s_data_v=std::shared_ptr<data_v >;
  using	datapair = std::pair<s_data_v,s_data_v>;
  const bool isGeno;
  size_t offset;
  std::unordered_map<T,s_data_v> data_buffer;
public:
  Genotype(const Rcpp::List data);
  Genotype(const Rcpp::List data,const Rcpp::List target,const bool isGenotype=true);
  void initialize(const std::pair<std::vector<T> ,std::vector<T> > &inp);
  datapair get(std::pair<T,T> index) const;
private:
  size_t validate(const Eigen::Map<Eigen::MatrixXd>  ref_data,const Eigen::Map<Eigen::ArrayXi> snp_index);
};





//#include "bits/input_bits.hpp"
