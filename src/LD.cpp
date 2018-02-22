#include <range/v3/all.hpp>
#include <LDshrink.h>
#include<algorithm>
#include <functional>
#include <tuple>
// [[Rcpp::depends(BH)]]
//[[Rcpp::plugins(cpp14)]]
#ifdef USE_MKL
#include "mkl.h"
#endif




//' Calculate the constant theta, given `m`
//' @param m a number indicating the size of the panel used to create the genetic map
//' (if using `1000-genomes-genetic-maps` from europeans, this number is 85)
//[[Rcpp::export(name="calc_theta")]]
double calc_theta_exp(const double m){
  return(calc_theta(m));
}



//' 'Melt' an LD matrix, dropping elements below a given r-square cutoff
//' @param m a number indicating the size of the panel used to create the genetic map
//' (if using `1000-genomes-genetic-maps` from europeans, this number is 85)
//[[Rcpp::export]]
Rcpp::DataFrame ld2df(const Matrix_external ldmat, Rcpp::StringVector rsid,const double r2cutoff=0.01){
  using namespace Rcpp;
  size_t p=ldmat.rows();
  if(p!=ldmat.cols()){
    Rcpp::stop("ldmat is not square!");
  }
  if(p!=rsid.size()){
    Rcpp::stop("rsid must be of length p!");
  }
  size_t dfsize = (p*(p-1))/2;
  std::vector<double>corv;
  corv.reserve(dfsize);
  std::vector<std::string> rowsnp;
  std::vector<std::string> colsnp;
  rowsnp.reserve(dfsize);
  colsnp.reserve(dfsize);

  // Eigen::ArrayXd corv(dfsize);
  // Rcpp::StringVector rowsnp(dfsize);
  // Rcpp::StringVector colsnp(dfsize);
  // Rcpp::Rcout<<"Generating DataFrame"<<std::endl;
  size_t k=0;
  for(int i=0; i<p;i++){
    for(int j=i+1; j<p;j++ ){
      double r2=ldmat.coeff(i,j)*ldmat.coeff(i,j);
      if(r2>r2cutoff){
        corv.push_back(r2);
        rowsnp.push_back(as<std::string>(rsid[i]));
        colsnp.push_back(as<std::string>(rsid[j]));
      }
    }
  }
  // Rcpp::Rcout<<"Returning DataFrame"<<std::endl;
  return(DataFrame::create(_["rowsnp"]=wrap(rowsnp),_["colsnp"]=wrap(colsnp),_["r2"]=wrap(corv),_["stringsAsFactors"]=false));
}


#ifdef USE_MKL
void calc_cov_mkl(Eigen::MatrixXd &X,Eigen::MatrixXd &S,Eigen::ArrayXd &meanv){
  VSLSSTaskPtr task;
  MKL_INT dim;
  MKL_INT n;
  MKL_INT x_storage;
  MKL_INT cov_storage;
  MKL_INT cor_storage;
  double *mean=meanv.data();
  double *x=X.data();
  double* cov=S.data();
  double* cor=nullptr;
  // cor[DIM][DIM];
  // double mean[DIM];
  double aTmp, rTmp;
  int i, j, errcode;
  int errnums = 0;
  //double pval_mean[DIM];
  //double pval_cov[DIM][DIM];
  /***** Initializing parameters for Summary Statistics task *****/
  dim         = X.cols();
  n           = X.rows();
  x_storage   = VSL_SS_MATRIX_STORAGE_COLS;
  cov_storage = VSL_SS_MATRIX_STORAGE_FULL;
  // for(i = 0; i < dim; i++)
  // {
  //   mean[i] = 0.0;
  //   for(j = 0; j < dim; j++)
  //   {
  //     cov[i][j] = 0;
  //     cor[i][j] = 0;
  //   }
  // }
  /***** Generate data set using VSL GaussianMV RNG *****/
 // errcode = dGenerateGaussianMVData( (double*)x, dim, n, (double*)a, (double*)C );
//CheckVslError(errcode);
  /***** Create Summary Statistics task *****/
  errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0, 0 );
  //CheckVslError(errcode);
  /***** Initialization of the task parameters using FULL_STORAGE
   for covariance/correlation matrices *****/
  errcode = vsldSSEditCovCor( task, mean, (double*)cov, &cov_storage,
                              NULL, NULL );
  //CheckVslError(errcode);
  /***** Compute covariance/correlation matrices using FAST method  *****/
  errcode = vsldSSCompute( task,
                           VSL_SS_COV,
                           VSL_SS_METHOD_FAST );
  errcode = vslSSDeleteTask( &task );
  //CheckVslError(errcode);
}
#endif

//[[Rcpp::export]]
Rcpp::NumericMatrix cov_mkl(Eigen::MatrixXd &X){
  const size_t p=X.cols();
  Eigen::MatrixXd S(p,p);
  Eigen::ArrayXd mean(p);
#ifdef USE_MKL
  calc_cov_mkl(X,S,mean);
  using namespace Rcpp;
  return(Rcpp::wrap(S));
#else
  Rcpp::stop("MKL not found!");
#endif
}



//[[Rcpp::export(name="calc_cov")]]
Eigen::MatrixXd calc_cov_exp(Eigen::MatrixXd &mat){
  Eigen::MatrixXd S;
  calc_cov(mat,S);
  return(S);
}


//[[Rcpp::export(name="cov_2_cor")]]
Eigen::MatrixXd cov_2_cor_exp(Eigen::MatrixXd &covmat){
  Eigen::MatrixXd tcovmat=covmat;
  cov_2_cor(tcovmat);
  return(tcovmat);
}


//[[Rcpp::export]]
Eigen::MatrixXd calcLD_prel(Eigen::MatrixXd hmata,const std::vector<double> mapa,const double m, const double Ne, const double cutoff){
  const size_t p=hmata.cols();
  Eigen::MatrixXd S(p,p);
  calcLD_pa(hmata,mapa,S,m,Ne,cutoff);
  return(S);
}


//[[Rcpp::export]]
Eigen::MatrixXd calcLD(Eigen::MatrixXd &hmata,const Eigen::ArrayXd &mapa,const double m, const double Ne, const double cutoff){

  double dosage_max=hmata.maxCoeff();
  bool isGeno=dosage_max>1;
  if(isGeno){
    Rcpp::Rcerr<<"LD is being estimated from genotype instead of haplotype"<<std::endl;
  }

  double theta=calc_theta(m);
  Eigen::MatrixXd S;
  calc_cov(hmata,S);
  if(isGeno){
    S*=0.5;
  }
  int numSNP=S.rows();
  S.triangularView<Eigen::StrictlyLower>().setZero();
  std::vector<double> mapv;
  mapv.resize(mapa.size());
  Eigen::VectorXd::Map(&mapv[0],mapa.size())=mapa;

  if(!is_sorted(mapv.begin(),mapv.end(),std::less<double>())){
    Rcpp::stop("Recombination map must be non-decreasing\n");
  }
  double tj=0;
  double ti=0;
  double rho=0;
  double tshrinkage;

  for(int i=0; i<numSNP;i++){
    ti=mapa(i);
    for(int j=i+1; j<numSNP;j++){
      tj=mapa(j);
      rho = 4*Ne*(tj-ti)/100;
      rho=-rho/(2*m);
      tshrinkage=std::exp(rho);
      if(tshrinkage<cutoff){
        tshrinkage=0;
      }
      S(i,j)=tshrinkage*S(i,j);
    }
  }
  // S=S+S.triangularView<Eigen::StrictlyUpper>().transpose();

  S.triangularView<Eigen::StrictlyLower>()=S.transpose();
  Eigen::ArrayXd eye(numSNP);
  eye.setOnes();
  eye=eye.array()*(0.5*theta * (1-0.5*theta));

  Eigen::MatrixXd SigHat = ((1-theta)*(1-theta))*S.array();
  SigHat.diagonal() = SigHat.diagonal().array()+eye;
  //       % SigHat is derived from Li and Stephens model (2003)
  //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  cov_2_cor(SigHat);
  return(SigHat);
}


//[[Rcpp::export]]
Rcpp::List calcLDt(Matrix_external hmata,arrayxd_external mapa,double m=85,double Ne=11490.672741,double cutoff=1e-3){
  using namespace Rcpp;
  double dosage_max=hmata.maxCoeff();
  bool isGeno=dosage_max>1;
  if(isGeno){
    Rcpp::Rcerr<<"LD is being estimated from genotype instead of haplotype"<<std::endl;
  }

  double theta=calc_theta(m);
  Eigen::MatrixXd S;
  calc_cov(hmata,S);
  if(isGeno){
    S*=0.5;
  }
  int numSNP=S.rows();
  S.triangularView<Eigen::StrictlyLower>().setZero();
  std::vector<double> mapv;
  mapv.resize(mapa.size());
  Eigen::VectorXd::Map(&mapv[0],mapa.size())=mapa;

  if(!is_sorted(mapv.begin(),mapv.end(),std::less<double>())){
    Rcpp::stop("Recombination map must be non-decreasing\n");
  }
  double tj=0;
  double ti=0;
  double rho=0;
  double tshrinkage;

  for(int i=0; i<numSNP;i++){
    ti=mapa(i);
    for(int j=i+1; j<numSNP;j++){
      tj=mapa(j);
      rho = 4*Ne*(tj-ti)/100;
      rho=-rho/(2*m);
      tshrinkage=std::exp(rho);
      if(tshrinkage<cutoff){
        tshrinkage=0;
      }
      S(i,j)=tshrinkage*S(i,j);
    }
  }

  // S=S+S.triangularView<Eigen::StrictlyUpper>().transpose();

  S.triangularView<Eigen::StrictlyLower>()=S.transpose();
  Eigen::ArrayXd eye(numSNP);
  eye.setOnes();
  eye=eye.array()*(0.5*theta * (1-0.5*theta));

  Eigen::MatrixXd SigHat = ((1-theta)*(1-theta))*S.array();
  SigHat.diagonal() = SigHat.diagonal().array()+eye;
  //       % SigHat is derived from Li and Stephens model (2003)
  //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  Eigen::MatrixXd tSigHat=SigHat;
  cov_2_cor(SigHat);
  return(Rcpp::List::create(_["SigHat"]=wrap(tSigHat),
                            _["R"]=wrap(SigHat),
                            _["S"]=wrap(S)));
}
