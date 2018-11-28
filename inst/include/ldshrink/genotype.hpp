#pragma once
#include "misc.hpp"
#include <RcppEigen.h>
#include <memory>
#include <Rinternals.h>

class Genotype{
  const int RTYPE= INTSXP;
  using T=int;
  using data_v = const std::pair<const Eigen::VectorXd,double>;
  //  using	s_data_v=data_v;
  using	datapair = std::pair<const data_v*,const data_v*>;
  const bool isGeno;
  size_t offset;
  std::unordered_map<T,data_v> data_buffer;
public:
  Genotype(const Rcpp::NumericMatrix data,const bool isGenotype=true);
  Genotype(const Rcpp::List data,const bool isGenotype=true);
  Genotype(const Rcpp::List data,const Rcpp::List target,const bool isGenotype=true);
  void initialize(const std::pair<std::vector<T> ,std::vector<T> > &inp);
  datapair get(std::array<T,2> index) const;
private:
  size_t validate(const Eigen::Map<Eigen::MatrixXd>  ref_data,const Eigen::Map<Eigen::ArrayXi> snp_index);
};





//#include "bits/input_bits.hpp"
