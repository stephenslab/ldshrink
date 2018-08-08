#include "LDshrink.h"
#include <algorithm>
#include <functional>
#include <tuple>

//[[Rcpp::depends(BH)]]

//[[Rcpp::export]]
SEXP sparse_LDshrink(Eigen::MatrixXd data, std::vector<double> mapd,Rcpp::IntegerVector indices,const double m, const double Ne, const double cutoff,const int total_size,const bool progress=true,const bool useLDshrink=true){

  const size_t p=mapd.size();
  std::vector<size_t> tindices(p);
  std::copy(indices.begin(),indices.end(),tindices.begin());
  LDshrink_data<double,Eigen::ColMajor> dat(std::move(data),std::move(mapd),std::move(tindices));
  dat.calc_vars();
  Sparse_cov<double,Eigen::ColMajor> sp(m,Ne,cutoff,0,progress);
  if(useLDshrink){
    sp.LDshrink(&dat);
  }else{
    sp.cor(&dat);
  }
  return(sp.dsCMatrix(total_size,total_size));
}


//[[Rcpp::export]]
SEXP sparse_LDshrink_p(Eigen::MatrixXd data_a,Eigen::MatrixXd data_b, std::vector<double> mapd_a,std::vector<double> mapd_b, Rcpp::IntegerVector indices_a,Rcpp::IntegerVector indices_b,const double m, const double Ne, const double cutoff,const int total_size,const bool progress=false,const bool useLDshrink = true){

  const size_t p_a=mapd_a.size();
  std::vector<size_t> tindices_a(p_a);
  std::copy(indices_a.begin(),indices_a.end(),tindices_a.begin());
  LDshrink_data<double,Eigen::ColMajor> dat_a(std::move(data_a),std::move(mapd_a),std::move(tindices_a));
  dat_a.calc_vars();
  const size_t p_b=mapd_b.size();
  std::vector<size_t> tindices_b(p_b);
  std::copy(indices_b.begin(),indices_b.end(),tindices_b.begin());
  LDshrink_data<double,Eigen::ColMajor> dat_b(std::move(data_b),std::move(mapd_b),std::move(tindices_b));
  dat_b.calc_vars();
  Sparse_cov<double,Eigen::ColMajor> sp(m,Ne,cutoff,0,progress);
  if(useLDshrink){
    sp.LDshrink(&dat_a,&dat_b);
  }else{
    sp.cor(&dat_a,&dat_b);
  }
  return(sp.dsCMatrix(total_size,total_size));
}


//[[Rcpp::export]]
Rcpp::DataFrame ld2df(Eigen::MatrixXd data, std::vector<double> mapd,Rcpp::RObject rsid,const double m, const double Ne, const double cutoff,const double r2cutoff=0.01,const bool progress=false,const bool useLDshrink = true){
  const size_t p=mapd.size();
  std::vector<size_t> indices(p);
  std::iota(indices.begin(), indices.end(), 0);
  if((rsid.sexp_type() != INTSXP) && (rsid.sexp_type() != STRSXP)){
    Rcpp::stop("rsid vector must be integer or string vector");
  }
  LDshrink_data<double,Eigen::ColMajor> dat(std::move(data),std::move(mapd),std::move(indices));
  dat.calc_vars();
  Sparse_cov<double,Eigen::ColMajor> sp(m,Ne,cutoff,r2cutoff,progress);
  if(useLDshrink){
    sp.LDshrink(&dat);
  }else{
    sp.cor(&dat);
  }
  if(rsid.sexp_type()==INTSXP){
    return(sp.ld2df(Rcpp::IntegerVector(rsid),Rcpp::IntegerVector(rsid)));
  }else{
    return(sp.ld2df(Rcpp::StringVector(rsid),Rcpp::StringVector(rsid)));
  }
}



//[[Rcpp::export]]
Rcpp::DataFrame ld2df_p(Eigen::MatrixXd data_a,Eigen::MatrixXd data_b, std::vector<double> mapd_a,std::vector<double> mapd_b, Rcpp::RObject rsid_a,Rcpp::RObject rsid_b,const double m, const double Ne, const double cutoff,const double r2cutoff=0.01,const bool progress=false,const bool useLDshrink = true){

  if(rsid_a.sexp_type()!=rsid_b.sexp_type()){
    Rcpp::stop("rsid_a and rsid_b must both be vectors of the same type");
  }
  if((rsid_a.sexp_type() != INTSXP) && (rsid_a.sexp_type() != STRSXP)){
    Rcpp::stop("rsid vectors must be integer or	string vectors (of the same type)");
  }

  const size_t p_a=mapd_a.size();
  std::vector<size_t> indices_a(p_a);
  std::iota(indices_a.begin(), indices_a.end(), 0);
  LDshrink_data<double,Eigen::ColMajor> dat_a(std::move(data_a),std::move(mapd_a),std::move(indices_a));
  dat_a.calc_vars();
  const size_t p_b=mapd_b.size();
  std::vector<size_t> indices_b(p_b);
  std::iota(indices_b.begin(), indices_b.end(), 0);
  LDshrink_data<double,Eigen::ColMajor> dat_b(std::move(data_b),std::move(mapd_b),std::move(indices_b));
  dat_b.calc_vars();
  Sparse_cov<double,Eigen::ColMajor> sp(m,Ne,cutoff,r2cutoff,progress);
  if(useLDshrink){
    sp.LDshrink(&dat_a,&dat_b);
  }else{
    sp.cor(&dat_a,&dat_b);
  }
  if(rsid_a.sexp_type()==INTSXP){
    return(sp.ld2df(Rcpp::IntegerVector(rsid_a),Rcpp::IntegerVector(rsid_b)));
  }else{
    return(sp.ld2df(Rcpp::StringVector(rsid_a),Rcpp::StringVector(rsid_b)));
  }
}
