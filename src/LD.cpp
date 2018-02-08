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
#include "boost/iterator/zip_iterator.hpp"


//

double calc_nmsum(const double m){
  int msize=(2*(int)m-1);
  Eigen::ArrayXd tx(msize);
  tx.setLinSpaced(msize,1,(int)(2*m-1));
  return  (1/tx).sum();
}


//[[Rcpp::export]]
bool less_pos(const int l_chrom,const int l_pos, const int r_chrom,const int r_pos){
  return(std::make_tuple(l_chrom,l_pos)<=std::make_tuple(r_chrom,r_pos));
}


struct SNPpos{
  std::pair<int,int> c_p;
  SNPpos(const int chrom,const int pos):c_p(std::make_pair(chrom,pos)){};
};


struct GRange{
  std::pair<SNPpos,SNPpos> range;
  GRange(const int chrom,const int start,const int end):range(std::make_pair(SNPpos(chrom,start),SNPpos(chrom,end))){};
};
bool operator<(const SNPpos p, const GRange g){
  return(p.c_p<g.range.first.c_p);
}
bool operator<=(const SNPpos p, const GRange g){
  return(p.c_p<=g.range.first.c_p);
}
bool operator>(const SNPpos p, const GRange g){
  return(p.c_p>g.range.second.c_p);
}
bool operator>=(const SNPpos p, const GRange g){
  return(p.c_p>=g.range.second.c_p);
}
bool operator==(const SNPpos p, const GRange g){
  return(p.c_p>=g.range.first.c_p & p.c_p<=g.range.second.c_p);
}

//[[Rcpp::export]]
bool sorted_snp_df(const Rcpp::DataFrame &snp_info){
  using namespace ranges;
  const Rcpp::IntegerVector chr= snp_info["chr"];
  const Rcpp::IntegerVector pos= snp_info["pos"];
  auto chr_r= make_iterator_range(chr.begin(),chr.end());
  auto pos_r= make_iterator_range(pos.begin(),pos.end());
  auto r = view::zip(chr_r,pos_r);
  bool is_sorted_res = ranges::is_sorted(r);
  return(is_sorted_res);
}


  


//[[Rcpp::export]]
Rcpp::IntegerVector set_ld_region(const Rcpp::DataFrame &ld_regions, const Rcpp::DataFrame &snp_info,const bool assign_all=true){

  const Rcpp::IntegerVector ld_chr= ld_regions["chr"];
  const Rcpp::IntegerVector ld_region_id= ld_regions["region_id"];
  const Rcpp::IntegerVector ld_start= ld_regions["start"];
  const Rcpp::IntegerVector ld_stop= ld_regions["stop"];
  const size_t ld_size=ld_chr.size();
  
  if(!std::is_sorted(ld_chr.begin(),ld_chr.end())){
    Rcpp::stop("break regions must be sorted by chromosome!");
  }

  const Rcpp::IntegerVector chr= snp_info["chr"];
  const Rcpp::IntegerVector pos= snp_info["pos"];
  const size_t snp_size = chr.size();
  Rcpp::IntegerVector ret_region(snp_size);
  if(!std::is_sorted(chr.begin(),chr.end())){
    Rcpp::stop("SNPs must be sorted by chromosome!");
  }
  size_t idx1=0;
  size_t idx2=0;
  
  while(idx1<snp_size){
    GRange gr(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
    //    auto range_start= std::make_tuple(ld_chr[idx2],ld_start[idx2]);
    //    auto range_stop= std::make_tuple(ld_chr[idx2],ld_stop[idx2]);
    const SNPpos s_pos(chr[idx1],pos[idx1]);
    if(s_pos==gr){
      //point in interval
      ret_region[idx1]=ld_region_id[idx2];
      idx1++;
    }else{
      if(s_pos<gr){
	// point before interval
        if(assign_all){
	  Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
	  Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
	}else{
	  ret_region[idx1]=NA_INTEGER;
	}
	idx1++;
      }
      if(s_pos>gr){
	//point after interval
	if(idx2<ld_size){
	  idx2++;
	}else{
	  if(assign_all){
	    Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
	    Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
	  }else{
	    ret_region[idx1]=NA_INTEGER;
	  }
          idx1++;
	}
      }
    }
  }
  return(ret_region);
}

//[[Rcpp::export]]
Rcpp::NumericVector interpolate_map(const Rcpp::NumericVector &map,const Rcpp::IntegerVector map_pos,const Rcpp::IntegerVector target_pos){
  
  const size_t map_snps=map.size();
  if(map_pos.size()!=map_snps){
    Rcpp::stop("map and map_pos must be the same size!");
  }
  const size_t t_snps=target_pos.size();
  Rcpp::NumericVector retvec(t_snps);
  
  size_t idx2=0;
  size_t idx1=0;
  while(idx1<t_snps){
    const int pos = target_pos[idx1];
    const int mappos = map_pos[idx2];
    if(pos == mappos){
      retvec[idx1]=map[idx2];
      idx1++;
    }else{
      if(pos<mappos){
        if(idx2==0){
          retvec[idx1]=map[idx2];
        }else{
          double prev_map = map[idx2-1];
          double prev_mappos = static_cast<double>(map_pos[idx2-1]);
          double frac = (static_cast<double>(pos)-prev_mappos)/(static_cast<double>(mappos)-prev_mappos);
          retvec[idx1]=prev_map+frac*(map[idx2]-prev_map);
          idx1++;
        }
      }else{
        if(pos>mappos){
          if(idx2==(map_snps-1)){
            retvec[idx1]=map[idx2];
            idx1++;
          }else{
            idx2++;
          }
        }
      }
    }
  }
  return(retvec);
}


//' Calculate the constant theta, given `m`
//' @param m a number indicating the size of the panel used to create the genetic map 
//' (if using `1000-genomes-genetic-maps` from europeans, this number is 85)
//[[Rcpp::export]]
double calc_theta(const double m){
  double nmsum=calc_nmsum(m);
  return((1/nmsum)/(2*m+1/nmsum));
}

void calc_cov( c_Matrix_internal mat,Eigen::MatrixXd &S){
  auto centered = mat.rowwise()-mat.colwise().mean();
  S= (((centered.adjoint()*centered)/double(mat.rows()-1)));  
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


void calc_cov_f(Eigen::MatrixXd &X){
  X=(((X.rowwise()-X.colwise().mean()).adjoint())*(X.rowwise()-X.colwise().mean()))/(X.rows()-1);
}

//[[Rcpp::export(name="calc_cov_f")]]
Eigen::MatrixXd calc_cov_exp_f(Eigen::MatrixXd mat){
  //Eigen::MatrixXd S;
  calc_cov_f(mat);
  return(mat);
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
  // 
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

//[[Rcpp::export]]
Rcpp::NumericMatrix cov_mkl(Eigen::MatrixXd &X){
  const size_t p=X.cols();
  Eigen::MatrixXd S(p,p);
  Eigen::ArrayXd mean(p);
  calc_cov_mkl(X,S,mean);
  using namespace Rcpp;
  return(Rcpp::wrap(S));
}

#endif




//[[Rcpp::export(name="calc_cov")]]
Eigen::MatrixXd calc_cov_exp(Matrix_external mat){
  Eigen::MatrixXd S;
  calc_cov(mat,S);
  return(S);
}


Eigen::ArrayXd calc_variance(c_Matrix_internal mat){
  int n=mat.rows();
  return(((mat.rowwise()-(mat.colwise().mean())).array().square().colwise().sum())/(n-1));
}

//[[Rcpp::export(name="calc_variance")]]
Eigen::ArrayXd calc_variance_exp(Matrix_external mat){
  return(calc_variance(mat));
}


void cov_2_cor(Matrix_internal covmat){
  Eigen::ArrayXd rowvar=1/covmat.diagonal().array().sqrt();
  covmat.array().colwise()*=rowvar;
  covmat.array().rowwise()*=rowvar.transpose();
  covmat.diagonal().setOnes();
}

//[[Rcpp::export(name="cov_2_cor")]]
Eigen::MatrixXd cov_2_cor_exp(Matrix_external covmat){
  Eigen::MatrixXd tcovmat=covmat;
  cov_2_cor(tcovmat);
  return(tcovmat);
  
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

  


void calcLD_pa(Eigen::MatrixXd &hmata,const std::vector<double> &mapa,Eigen::MatrixXd & S,const double m, const double Ne, const double cutoff){
  
  double dosage_max=hmata.maxCoeff();
  bool isGeno=dosage_max>1;
  if(isGeno){
    //Rcpp::Rcerr<<"LD is being estimated from genotype instead of haplotype"<<std::endl;
  }
  
  double theta=calc_theta(m);
  calc_cov(hmata,S);
  if(isGeno){
    S*=0.5;
  }
  int numSNP=S.rows();
  S.triangularView<Eigen::StrictlyLower>().setZero();

  
  if(!is_sorted(mapa.begin(),mapa.end(),std::less<double>())){
    
    Rcpp::stop("Recombination map must be non-decreasing\n");
  }
  double tj=0;
  double ti=0;
  double rho=0;
  double tshrinkage;
  
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
  
  S.triangularView<Eigen::StrictlyLower>()=S.transpose();
  Eigen::ArrayXd eye(numSNP);
  eye.setOnes();
  eye=eye.array()*(0.5*theta * (1-0.5*theta));
  
  S = ((1-theta)*(1-theta))*S.array();
  S.diagonal() = S.diagonal().array()+eye;
  //       % SigHat is derived from Li and Stephens model (2003)
  //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  cov_2_cor(S);
 // return(SigHat);
}

//[[Rcpp::export]]
Eigen::MatrixXd calcLD_prel(Eigen::MatrixXd hmata,const std::vector<double> mapa,const double m, const double Ne, const double cutoff){
  const size_t p=hmata.cols();
  Eigen::MatrixXd S(p,p);
  calcLD_pa(hmata,mapa,S,m,Ne,cutoff);
  return(S);
}

Eigen::MatrixXd calcLD(const c_Matrix_internal hmata,const c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff){
  
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
























//[[Rcpp::export(name="calcLD")]]
Eigen::MatrixXd calcLD_exp(Matrix_external hmata,arrayxd_external mapa,const double m=85, const double Ne=11490.672741, const double cutoff=1e-3){
  return(calcLD(hmata,mapa,m,Ne,cutoff));
}



Eigen::SparseMatrix<double> sp_calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff){
  
  // compute_shrinkage_cor(dist,S,hmata,theta,m,Ne,cutoff);
  
  Eigen::MatrixXd dist = calcLD(hmata,mapa,m,Ne,cutoff);
  // dist.triangularView<Eigen::Lower>()=dist.transpose();
  Eigen::SparseMatrix<double> retmat=dist.sparseView();
  return(retmat);
}

//[[Rcpp::export(name="sp_calcLD")]]
Eigen::SparseMatrix<double> sp_calcLD_exp(Matrix_external hmata,arrayxd_external mapa,const double m, const double Ne, const double cutoff){
  return(sp_calcLD(hmata,mapa,m,Ne,cutoff));
}



Eigen::SparseMatrix<double> sp_calcLD_symm(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff){
  
  // compute_shrinkage_cor(dist,S,hmata,theta,m,Ne,cutoff);
  
  Eigen::MatrixXd dist = calcLD(hmata,mapa,m,Ne,cutoff);
  // dist.triangularView<Eigen::Lower>()=dist.transpose();
  Eigen::SparseMatrix<double> retmat=Eigen::MatrixXd(dist.triangularView<Eigen::Upper>()).sparseView();
  return(retmat);
}

//[[Rcpp::export(name="sp_calcLD_symm")]]
Eigen::SparseMatrix<double> sp_calcLD_symm_exp(Matrix_external hmata,arrayxd_external mapa,const double m, const double Ne, const double cutoff){
  return(sp_calcLD_symm(hmata,mapa,m,Ne,cutoff));
}
