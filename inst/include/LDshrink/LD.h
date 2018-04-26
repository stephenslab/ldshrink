#ifndef LD_H
#define LD_H
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <memory>
#include "LDshrink_types.h"



class IO_h5{
  friend class MatSlices;
protected:
  std::unordered_map<std::string,std::shared_ptr<HighFive::File> >  m_file_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::Group> >  m_group_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::DataSet> > m_dataset_map;
public:
  MatSlices input_f;
  MatSlices output_f;
  IO_h5(const Rcpp::DataFrame input_dff,
       const Rcpp::DataFrame output_dff):m_file_map(),
       m_group_map(),
       m_dataset_map(),
       input_f(input_dff,m_file_map,m_group_map,m_dataset_map,true),
       output_f(output_dff,m_file_map,m_group_map,m_dataset_map,false){}
  size_t num_input_groups()const {
    return(input_f.chunk_map.size());
  }
  size_t num_output_groups()const{
    return(output_f.chunk_map.size());
  }
};


//class IO_R{



/* template<typename IOT> */
/* class LDshrink{ */




/*   /\* std::unordered_map<int,std::shared_future<mmat> > Hmap; *\/ */
/*   /\* std::unordered_map<int,std::shared_future<mmat> > Smap; *\/ */
/*   /\* std::unordered_map<int,std::shared_future<mmat> > Dmap; *\/ */
/*   /\* std::unordered_map<int,std::shared_future<mmat> > Qmap; *\/ */
/*   mmati iv; */
/*   IOT io_obj; */
/*   const size_t num_reg; */
/*   const double m; */
/*   const double Ne; */
/*   const double cutoff; */
/*   const bool SNPfirst; */
/*   const bool doEVD; */
/*   const bool doSVD; */
/*   const bool doDF; */
/*   const double r2cutoff; */


/*   mmat Rsq; */


/*  public: */
/*  LDshrink(IOT & io_obj_, */
/* 	  const double m_, */
/* 	  const double Ne_, */
/* 	  const double cutoff_, */
/* 	  const bool SNPfirst_=true, */
/* 	  const bool doEVD_=false, */
/* 	  const bool doSVD_=false, */
/* 	  const bool doDF_=false, */
/* 	  const double r2cutoff_=0.01):io_obj(io_obj_), */
/*     num_reg(io_obj.num_input_groups()), */
/*     m(m_), */
/*     Ne(Ne_), */
/*     cutoff(cutoff_), */
/*     SNPfirst(SNPfirst_), */
/*     doEVD(doEVD_), */
/*     doSVD(doSVD_), */
/*     doDF(doDF_), */
/*     r2cutoff(r2cutoff_), */
/*     prog_bar(num_reg,true){ */


/*   } */
/*   void calcLD(){ */



/*   } */

/* }; */




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



template<typename T>
class LDshrinker{
  std::vector<T> Sdat;
public:
  const T m;
  const T Ne;
  const T theta;
  const T cutoff;
  Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > S;
  Eigen::Map<const Eigen::Array<T,Eigen::Dynamic,1> > mapd;
  //  const size_t N;
  LDshrinker(const T m_,const T Ne_,const T cutoff_):m(m_),
						   theta(calc_theta(m)),
						   cutoff(cutoff_),
						     Ne(Ne_),
						     mapd(nullptr,1),
						     S(nullptr,1,1)
  {}
  LDshrinker(const T m_,const T Ne_,const T cutoff_,
	     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &hmat_,
	     Eigen::Map<const Eigen::Array<T,Eigen::Dynamic,1> > &map_):m(m_),
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
	     Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > &hmat_,
	     Eigen::Map<const Eigen::Array<T,Eigen::Dynamic,1> > &map_):m(m_),
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
  

  LDshrinker(Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > &S_,
	     Eigen::Map<const Eigen::Array<T,Eigen::Dynamic,1> >  &map_,const T m_,const T Ne_,const T cutoff_):m(m_),
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

  static auto calc_cov(  Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > &hmata){
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



#endif
