#include "ldshrink.hpp"
#include <Rcpp.h>

// The Genotype class acts as an intermediary between the data (either on disk or in R), and the ldshrink algorithm
// This means doing any	precomputations	necessary (including subsetting rows or columns if provided)
Genotype::Genotype(const Rcpp::List data,const bool isGenotype):isGeno(isGenotype),offset(0){
  using namespace Rcpp;
  if (data.size() == 2) {
      this->validate(as<Eigen::Map<Eigen::MatrixXd> >(data["data"]),as<Eigen::Map<Eigen::ArrayXi>>(data["snp_index"]));
  } else {
    this->validate(
        as<Eigen::Map<Eigen::MatrixXd>>(as<List>(data[0])["data"]),
        as<Eigen::Map<Eigen::ArrayXi>>(as<List>(data[0])["snp_index"]));
    this->validate(
        as<Eigen::Map<Eigen::MatrixXd>>(as<List>(data[1])["data"]),
        as<Eigen::Map<Eigen::ArrayXi>>(as<List>(data[1])["snp_index"]));
  }
}

Genotype::Genotype(const Rcpp::NumericMatrix data,const bool isGenotype):isGeno(isGenotype),offset(0){
  using namespace Rcpp;
  auto data_mat = as<Eigen::Map<Eigen::MatrixXd>>(data);
  const size_t max_p=data_mat.cols();
  Eigen::ArrayXi data_ind =Eigen::ArrayXi::LinSpaced(max_p,1,max_p);
  auto mapi= Eigen::Map<Eigen::ArrayXi>(data_ind.data(),data_ind.size());
  this->validate(data_mat, mapi);

}

Genotype::Genotype(const Rcpp::List data,const Rcpp::List target,const bool isGenotype):isGeno(isGenotype),offset(0){
  using namespace Rcpp;
  //  if (data.size() == 1) {
  auto data_mat =as<Eigen::Map<Eigen::MatrixXd>>(data["data"]);
  const size_t max_p=data_mat.cols();
  auto data_ind =value_or<Eigen::ArrayXi>(target,"snp_index",Eigen::ArrayXi::LinSpaced(max_p,1,max_p));
  auto mapi= Eigen::Map<Eigen::ArrayXi>(data_ind.data(),data_ind.size());
  this->validate(data_mat,mapi);
  // } else {
  //   this->validate(
  //       as<Eigen::Map<Eigen::MatrixXd>>(as<List>(data[1])["data"]),
  //       as<Eigen::Map<Eigen::ArrayXi>>(as<List>(target[1])["snp_index"]));
  //   this->validate(
  //       as<Eigen::Map<Eigen::MatrixXd>>(as<List>(data[0])["data"]),
  //       as<Eigen::Map<Eigen::ArrayXi>>(as<List>(target[0])["snp_index"]));
  // }
}

void Genotype::initialize(const std::pair<std::vector<T> ,std::vector<T> > &inp){

}


size_t Genotype::validate(const Eigen::Map<Eigen::MatrixXd>  ref_data,const Eigen::Map<Eigen::ArrayXi> snp_index){
  using namespace Rcpp;


  const size_t N= ref_data.rows();
  const size_t p= snp_index.size();
  const double genoflip	= isGeno ? 2 : 1;
  Eigen::VectorXd tempd(N);
  const size_t t_offset = offset;
  for(size_t i=0; i<p; i++){
    double sum_d=0;
    const int snpi= snp_index(i)-1;
    for (size_t j = 0; j < N; j++) {
      tempd = ref_data.col(snpi);
    }
    tempd = tempd.array() - tempd.mean();
    data_buffer.emplace(std::make_pair(
        i + t_offset,
        std::make_shared<data_v>(tempd, tempd.array().square().sum() /
                                            static_cast<double>(N - 1))));
    offset++;
  }
  return (offset);
}

typename Genotype::datapair Genotype::get(std::pair<T,T> index) const {
  using ret_datapair = typename Genotype::datapair;
  return(std::make_pair(data_buffer.find(index.first)->second,data_buffer.find(index.second)->second));
}
