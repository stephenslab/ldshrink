#ifndef LD_H
#define LD_H
#include <RcppEigen.h>


#include "LDshrink_types.h"


inline double calc_nmsum(const double m){
  int msize=(2*(int)m-1);
  Eigen::ArrayXd tx(msize);
  tx.setLinSpaced(msize,1,(int)(2*m-1));
  return  (1/tx).sum();
}

inline double calc_theta(const double m){
  double nmsum=calc_nmsum(m);
  return((1/nmsum)/(2*m+1/nmsum));
}

template<typename DerivedA,typename DerivedB>
  void calc_cov( Eigen::MatrixBase<DerivedA> &mat,Eigen::MatrixBase<DerivedB> &S){
   auto centered = mat.rowwise()-mat.colwise().mean();
  S= (((centered.adjoint()*centered)/(mat.rows()-1)));
}
template<typename DerivedA,typename DerivedB>
void calc_cov_s(Eigen::MatrixBase<DerivedA> &mat,Eigen::MatrixBase<DerivedB> &S){
  auto n = Eigen::nbThreads( );
  // Rcpp::Rcout<<"Using : "<<n<<" threads"<<std::endl;
  // mat = mat.rowwise()-mat.colwise().mean();
  S= (((mat.adjoint()*mat)/(mat.rows()-1)));
}


template<typename Derived,int Flags=Eigen::internal::traits<Derived>::Flags & Eigen::RowMajorBit ? Eigen::RowMajor : Eigen::ColMajor>
void cov_2_cor(Eigen::MatrixBase<Derived> &covmat){
  typedef typename Derived::Scalar Scalar;
  Eigen::Array<Scalar,Eigen::Dynamic,1> rowvar = 1/covmat.diagonal().array().sqrt();
  covmat.array().colwise()*=rowvar;
  covmat.array().rowwise()*=rowvar.transpose();
  covmat.diagonal().setOnes();
}

template<typename Derived,int Flags=Eigen::internal::traits<Derived>::Flags & Eigen::RowMajorBit ? Eigen::RowMajor : Eigen::ColMajor>
void calcLD_pa(Eigen::MatrixBase<Derived> &hmata,
               const std::vector<typename Derived::Scalar> &mapa,
               Eigen::MatrixBase<Derived> &S,
               const  typename Derived::Scalar m,
               const typename Derived::Scalar Ne,
               const typename Derived::Scalar cutoff){
  typedef typename Derived::Scalar Scalar;
  Scalar dosage_max=hmata.maxCoeff();
  bool isGeno=dosage_max>1;
  if(isGeno){
    //Rcpp::Rcerr<<"LD is being estimated from genotype instead of haplotype"<<std::endl;
  }
  Scalar theta=calc_theta(m);
  hmata = hmata.rowwise()-hmata.colwise().mean();
  calc_cov_s(hmata,S);
  if(isGeno){
    S*=0.5;
  }
  int numSNP=S.rows();
  auto nts=S.template triangularView<Eigen::StrictlyLower>().setZero();

  if(!is_sorted(mapa.begin(),mapa.end(),std::less<Scalar>())){
    Rcpp::stop("Recombination map must be non-decreasing\n");
  }
  Scalar tj=0;
  Scalar ti=0;
  Scalar rho=0;
  Scalar tshrinkage;

  for(int i=0; i<numSNP;i++){
    ti=mapa[i];
    for(int j=i+1; j<numSNP;j++){
      tj=mapa[j];
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
  S.template triangularView<Eigen::StrictlyLower>()=S.transpose();
  // Eigen::Matrix<Scalar,Eigen::Dynamic,1,Flags> eyem(numSNP);

  S = ((1-theta)*(1-theta))*S.array();
  for(int i=0; i<numSNP;i++){
    S(i,i)+=0.5*theta * (1-0.5*theta);
  }
  // S.diagonal()+=eyem.asDiagonal();
  //       % SigHat is derived from Li and Stephens model (2003)
  //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  cov_2_cor(S);
 // return(SigHat);
}

template<typename Derived>
void LD2df(const Eigen::MatrixBase<Derived> &ldmat,
           const std::vector<std::string> &rsid,
           std::vector<std::string> &rowsnp,
           std::vector<std::string> &colsnp,
           std::vector<double> &corv,
           const double r2cutoff=0.01){
  size_t p=ldmat.rows();
  if(p!=ldmat.cols()){
    Rcpp::stop("ldmat is not square!");
  }
  if(p!=rsid.size()){
    Rcpp::stop("rsid must be of length p!");
  }
  size_t dfsize = (p*(p-1))/2;
  if(corv.capacity()<dfsize){
    corv.reserve(dfsize);
  }
  if(rowsnp.capacity()<dfsize){
    rowsnp.reserve(dfsize);
  }
  if(colsnp.capacity()<dfsize){
    colsnp.reserve(dfsize);
  }
  
  
  size_t k=0;
  for(int i=0; i<p;i++){
    for(int j=i+1; j<p;j++ ){
      const double r2=ldmat.coeff(i,j)*ldmat.coeff(i,j);
      if(r2>r2cutoff){
        corv.push_back(r2);
        rowsnp.push_back(rsid[i]);
        colsnp.push_back(rsid[j]);
      }
    }
  }
}



#endif
