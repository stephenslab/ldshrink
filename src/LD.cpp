#include <LDshrink.h>
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
  Eigen::Map<const Eigen::ArrayXd> mmap(&mapd(0),p);
  LDshrinker<double> lds(m,Ne,cutoff,mH,mmap);
  lds.Shrink();
  return(Rcpp::wrap(lds.S));
}

//[[Rcpp::export]]
Rcpp::NumericMatrix shrinkCov(const Rcpp::NumericMatrix S,const Rcpp::NumericVector &mapd,const double m, const double Ne, const double cutoff){
  const size_t p = mapd.size();
  Rcpp::NumericMatrix nS(Rcpp::clone(S));
  Eigen::Map<Eigen::MatrixXd> mS(&nS(0,0),p,p);
  Eigen::Map<const Eigen::Array<double,Eigen::Dynamic,1> > mmap(&mapd(0),p);
  Eigen::Map<const Eigen::Array<double,Eigen::Dynamic,1> > tmmap(mmap.data(),p);

  LDshrinker<double> lds(mS,mmap,m,Ne,cutoff);
  lds.Shrink();
  return(nS);
}



//[[Rcpp::export]]
Eigen::MatrixXd calcDist(Eigen::ArrayXd &map){
  const size_t p=map.size();
  return(  ((map.transpose().colwise().replicate(p)).colwise()-map));
}





// Eigen::MatrixXd calcLD(Eigen::MatrixXd &hmata,const Eigen::ArrayXd &mapa,const double m, const double Ne, const double cutoff){
//   double dosage_max=hmata.maxCoeff();
//   bool isGeno=dosage_max>1;
//   if(isGeno){
//     Rcpp::Rcerr<<"LD is being estimated from genotype instead of haplotype"<<std::endl;
//   }

//   double theta=calc_theta(m);
//   Eigen::MatrixXd S;
//   calc_cov(hmata,S);
//   if(isGeno){
//     S*=0.5;
//   }
//   int numSNP=S.rows();
//   S.triangularView<Eigen::StrictlyLower>().setZero();
//   std::vector<double> mapv;
//   mapv.resize(mapa.size());
//   Eigen::VectorXd::Map(&mapv[0],mapa.size())=mapa;

//   if(!is_sorted(mapv.begin(),mapv.end(),std::less<double>())){
//     Rcpp::stop("Recombination map must be non-decreasing\n");
//   }
//   double tj=0;
//   double ti=0;
//   double rho=0;
//   double tshrinkage;

//   for(int i=0; i<numSNP;i++){
//     ti=mapa(i);
//     for(int j=i+1; j<numSNP;j++){
//       tj=mapa(j);
//       rho = 4*Ne*(tj-ti)/100;
//       rho=-rho/(2*m);
//       tshrinkage=std::exp(rho);
//       if(tshrinkage<cutoff){
//         tshrinkage=0;
//       }
//       S(i,j)=tshrinkage*S(i,j);
//     }
//   }
//   // S=S+S.triangularView<Eigen::StrictlyUpper>().transpose();

//   S.triangularView<Eigen::StrictlyLower>()=S.transpose();
//   Eigen::ArrayXd eye(numSNP);
//   eye.setOnes();
//   eye=eye.array()*(0.5*theta * (1-0.5*theta));

//   Eigen::MatrixXd SigHat = ((1-theta)*(1-theta))*S.array();
//   SigHat.diagonal() = SigHat.diagonal().array()+eye;
//   //       % SigHat is derived from Li and Stephens model (2003)
//   //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
//   // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
//   cov_2_cor(SigHat);
//   return(SigHat);
// }



// Rcpp::List calcLDt(Matrix_external hmata,arrayxd_external mapa,double m=85,double Ne=11490.672741,double cutoff=1e-3){
//   using namespace Rcpp;
//   double dosage_max=hmata.maxCoeff();
//   bool isGeno=dosage_max>1;
//   if(isGeno){
//     Rcpp::Rcerr<<"LD is being estimated from genotype instead of haplotype"<<std::endl;
//   }

//   double theta=calc_theta(m);
//   Eigen::MatrixXd S;
//   calc_cov(hmata,S);
//   if(isGeno){
//     S*=0.5;
//   }
//   int numSNP=S.rows();
//   S.triangularView<Eigen::StrictlyLower>().setZero();
//   std::vector<double> mapv;
//   mapv.resize(mapa.size());
//   Eigen::VectorXd::Map(&mapv[0],mapa.size())=mapa;

//   if(!is_sorted(mapv.begin(),mapv.end(),std::less<double>())){
//     Rcpp::stop("Recombination map must be non-decreasing\n");
//   }
//   double tj=0;
//   double ti=0;
//   double rho=0;
//   double tshrinkage;

//   for(int i=0; i<numSNP;i++){
//     ti=mapa(i);
//     for(int j=i+1; j<numSNP;j++){
//       tj=mapa(j);
//       rho = 4*Ne*(tj-ti)/100;
//       rho=-rho/(2*m);
//       tshrinkage=std::exp(rho);
//       if(tshrinkage<cutoff){
//         tshrinkage=0;
//       }
//       S(i,j)=tshrinkage*S(i,j);
//     }
//   }

//   // S=S+S.triangularView<Eigen::StrictlyUpper>().transpose();

//   S.triangularView<Eigen::StrictlyLower>()=S.transpose();
//   Eigen::ArrayXd eye(numSNP);
//   eye.setOnes();
//   eye=eye.array()*(0.5*theta * (1-0.5*theta));

//   Eigen::MatrixXd SigHat = ((1-theta)*(1-theta))*S.array();
//   SigHat.diagonal() = SigHat.diagonal().array()+eye;
//   //       % SigHat is derived from Li and Stephens model (2003)
//   //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
//   // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
//   Eigen::MatrixXd tSigHat=SigHat;
//   cov_2_cor(SigHat);
//   return(Rcpp::List::create(_["SigHat"]=wrap(tSigHat),
//                             _["R"]=wrap(SigHat),
//                             _["S"]=wrap(S)));
// }
