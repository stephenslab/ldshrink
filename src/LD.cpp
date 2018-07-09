#include "LDshrink.h"
#include<algorithm>
#include <functional>
#include <tuple>
// [[Rcpp::depends(BH)]]


//' Calculate the constant theta, given `m`
//' @param m a number indicating the size of the panel used to create the genetic map
//' (if using `1000-genomes-genetic-maps` from europeans, this number is 85)
//[[Rcpp::export(name="calc_theta")]]
double calc_theta_exp(const double m){
  return(calc_theta(m));
}






//' 'Melt' an LD matrix, dropping elements below a given r-square cutoff
//[[Rcpp::export]]
Rcpp::DataFrame ld2df(const Eigen::MatrixXd &ldmat, Rcpp::StringVector rsid,const double r2cutoff=0.01,const bool stringsAsFactors=false){
  using namespace Rcpp;
  std::vector<std::string> rowsnp;
  std::vector<std::string> colsnp;
  std::vector<double>corv;
  LD2df(ldmat,Rcpp::as<std::vector<std::string> >(rsid),rowsnp,colsnp,corv,r2cutoff);
  return(Rcpp::DataFrame::create(_["rowsnp"]=Rcpp::wrap(rowsnp),_["colsnp"]=Rcpp::wrap(colsnp),_["r2"]=wrap(corv),_["stringsAsFactors"]=stringsAsFactors));
}


//[[Rcpp::export]]
Rcpp::NumericMatrix shrinkCov(const Rcpp::NumericMatrix S,const Rcpp::NumericVector &mapd,const double m, const double Ne, const double cutoff){
  const size_t p = mapd.size();
  Rcpp::NumericMatrix nS(Rcpp::clone(S));
  Eigen::Map<Eigen::MatrixXd> mS(&nS(0,0),p,p);
  Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > mmap(&mapd(0),p);
  // Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > tmmap(mmap.data(),p);
  LDshrinker<double> lds(mS,mmap,m,Ne,cutoff);
  lds.Shrink();
  return(nS);
}

//[[Rcpp::export]]
Rcpp::NumericMatrix fastLDshrink(const Rcpp::NumericMatrix genotype_data,const Rcpp::NumericVector &mapd,const double m, const double Ne, const double cutoff,const bool isGeno=true, const bool cov_2_cor=true){
  const size_t p = mapd.size();
  Rcpp::NumericMatrix nS(p,p);
  Eigen::Map<Eigen::MatrixXd> mS(&nS(0,0),p,p);
  Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > mmap(&mapd(0),p);
  LDshrinker<double> lds(mS,mmap,m,Ne,cutoff,genotype_data,isGeno);
  // Rcpp::Rcerr<<Eigen::nbThreads( )<<std::endl;
  lds.Shrink();
  if(cov_2_cor){
    lds.cov_2_cor();
  }
  return(nS);
}





//[[Rcpp::export]]
Eigen::MatrixXd calcDist(Eigen::ArrayXd &map){
  const size_t p=map.size();
  return(  ((map.transpose().colwise().replicate(p)).colwise()-map));
}

