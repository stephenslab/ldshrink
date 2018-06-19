#pragma once
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <memory>
#include "mapdiff.hpp"




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
  void calc_cov( Eigen::MatrixBase<DerivedA> mat,Eigen::MatrixBase<DerivedB> &S){
    mat = mat.rowwise()-mat.colwise().mean();
  S= (((mat.adjoint()*mat)/(mat.rows()-1)));
}
template<typename DerivedA,typename DerivedB>
void calc_cov_s(Eigen::MatrixBase<DerivedA> &mat,Eigen::MatrixBase<DerivedB> &S){
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



template<typename T,int RM=Eigen::ColMajor>
class LDshrinker{
  std::vector<T> Sdat;
public:
  const T m;
  const T Ne;
  const T theta;
  const T cutoff;
  Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,RM> > S;
  Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1> > mapd;
  //  const size_t N;
  LDshrinker(const T m_,const T Ne_,const T cutoff_):m(m_),
						   theta(calc_theta(m)),
						   cutoff(cutoff_),
						     Ne(Ne_),
						     mapd(nullptr,1),
						     S(nullptr,1,1)
  {}
  LDshrinker(const T m_,const T Ne_,const T cutoff_,
	     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,RM> &hmat_,
	     Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1> > &map_):m(m_),
									theta(calc_theta(m)),
									cutoff(cutoff_),
									Ne(Ne_),
									mapd(map_.data(),map_.size()),
  									S(nullptr,map_.size(),map_.size()){
    if(!is_sorted(mapd.data(),mapd.data()+mapd.size(),std::less<T>())){
      Rcpp::stop("Recombination map must be non-decreasing\n");
    }
    const size_t p=mapd.size();
    if(hmat_.cols()!=p){
      Rcpp::stop("Reference panel data matrix must have column number equal to the size of genetic map");
    }
    Sdat.resize(p*p);
    S=Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> >(Sdat.data(),p,p);
    S=calc_cov(hmat_);
  }

  LDshrinker(const T m_,const T Ne_,const T cutoff_,
	     Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,RM> > &hmat_,
	     Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1> > &map_):m(m_),
									theta(calc_theta(m)),
									cutoff(cutoff_),
									Ne(Ne_),
									mapd(map_.data(),map_.size()),
									S(nullptr,map_.size(),map_.size())
  {
    if(!is_sorted(mapd.data(),mapd.data()+mapd.size(),std::less<T>())){
      Rcpp::stop("Recombination map must be non-decreasing\n");
    }
    const size_t p=mapd.size();
    if(hmat_.cols()!=p){
      Rcpp::stop("Reference panel data matrix must have column number equal to the size of genetic map");
    }
    Sdat.resize(p*p);
    S=Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> >(Sdat.data(),p,p);
    S=calc_cov(hmat_);
  }
  

  LDshrinker(Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,RM> > &S_,
	     Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1> >  &map_,const T m_,const T Ne_,const T cutoff_):m(m_),
										    theta(calc_theta(m)),
										    cutoff(cutoff_),
										    Ne(Ne_),
										    mapd(map_.data(),map_.size()),
										    S(S_.data(),S_.rows(),S_.cols())
  {
    if(!is_sorted(mapd.data(),mapd.data()+mapd.size(),std::less<T>())){
      Rcpp::stop("Recombination map must be non-decreasing\n");
    }
    const size_t p=mapd.size();
    if(S.cols()!=p || S.rows()!=p){
      Rcpp::stop("Covariance Matrix have column(row) number equal to the size of genetic map");
    }
  }

  static auto calc_cov(  Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,RM> > &hmata){
    T dosage_max=hmata.maxCoeff();
    bool isGeno=dosage_max>1;
    const double GenoMult = isGeno ? 0.5 : 1;
    hmata = hmata.rowwise()-hmata.colwise().mean();
    return((((hmata.adjoint()*hmata)/(hmata.rows()-1)))*GenoMult);
  }

  T shrinkp(const T map_dist){
    T rho = 4*Ne*(map_dist)/100;
    rho=-rho/(2*m);
    rho =std::exp(rho);
    return(rho<cutoff ? 0 : rho);
  }


  void Shrink(){
    //    const size_t p=map.size();
    const int numSNP=mapd.size();
    T ti=0;
    T tshrinkage;
    for(int i=0; i<numSNP;i++){
      ti=mapd[i];
      for(int j=i+1; j<numSNP;j++){
	S(i,j)=S(i,j)*shrinkp(mapd[j]-ti);
      }
    }

    S.template triangularView<Eigen::StrictlyLower>()=S.transpose();
    S = ((1-theta)*(1-theta))*S.array();
    for(int i=0; i<numSNP;i++){
      S(i,i)+=0.5*theta * (1-0.5*theta);
    }
  }


  // void alt_Shrink(){
  //   const int numSNP=mapd.size();
  //   S=S*makeMapDiff(mapd,m,Ne,cutoff);
  //   S = ((1-theta)*(1-theta))*S.array();
  //   for(int i=0; i<numSNP;i++){
  //     S(i,i)+=0.5*theta * (1-0.5*theta);
  //   }
  // }



};


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
