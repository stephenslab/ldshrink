#include "LDshrink.h"
#include<algorithm>
#include <functional>
#include <tuple>
// [[Rcpp::depends(BH)]]
//[[Rcpp::plugins(cpp14)]]


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
  return(DataFrame::create(_["rowsnp"]=wrap(rowsnp),_["colsnp"]=wrap(colsnp),_["r2"]=wrap(corv),_["stringsAsFactors"]=stringsAsFactors));
}

//[[Rcpp::export]]
Rcpp::NumericMatrix shrinkPanel(const Rcpp::NumericMatrix hpanel,const Rcpp::NumericVector &mapd,const double m, const double Ne, const double cutoff){
  const size_t p = mapd.size();
  const size_t n=hpanel.rows();
  if(hpanel.ncol()!=p){
    Rcpp::stop("Reference panel data matrix must have column number equal to the size of genetic map");
  }
  Rcpp::NumericMatrix nH(Rcpp::clone(hpanel));
  Eigen::Map<Eigen::MatrixXd> mH(&nH(0,0),n,p);
  Eigen::Map<const Eigen::VectorXd> mmap(&mapd(0),p);
  LDshrinker<double> lds(m,Ne,cutoff,mH,mmap);
  lds.Shrink();
  return(Rcpp::wrap(lds.S));
}

//[[Rcpp::export]]
Rcpp::NumericMatrix shrinkCov(const Rcpp::NumericMatrix S,const Rcpp::NumericVector &mapd,const double m, const double Ne, const double cutoff){
  const size_t p = mapd.size();
  Rcpp::NumericMatrix nS(Rcpp::clone(S));
  Eigen::Map<Eigen::MatrixXd> mS(&nS(0,0),p,p);
  Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > mmap(&mapd(0),p);
  Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > tmmap(mmap.data(),p);

  LDshrinker<double> lds(mS,mmap,m,Ne,cutoff);
  lds.Shrink();
  return(nS);
}



// Rcpp::NumericMatrix alt_shrinkCov(const Rcpp::NumericMatrix S,const Rcpp::NumericVector &mapd,const double m, const double Ne, const double cutoff){
//   const size_t p = mapd.size();
//   Rcpp::NumericMatrix nS(Rcpp::clone(S));
//   Eigen::Map<Eigen::MatrixXd> mS(&nS(0,0),p,p);
//   Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > mmap(&mapd(0),p);
//   Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1> > tmmap(mmap.data(),p);

//   LDshrinker<double> lds(mS,mmap,m,Ne,cutoff);
//   lds.alt_Shrink();
//   return(nS);
// }





//[[Rcpp::export]]
Eigen::MatrixXd calcDist(Eigen::ArrayXd &map){
  const size_t p=map.size();
  return(  ((map.transpose().colwise().replicate(p)).colwise()-map));
}
