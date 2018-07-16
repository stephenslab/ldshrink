#include "LDshrink.h"
#include <algorithm>
#include <functional>
#include <tuple>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <tbb/tbb.h>
//[[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


std::pair<size_t,size_t> row_col_t(const size_t k, const size_t p){
  size_t i = floor( ( 2*p+1 - sqrt( (2*p+1)*(2*p+1) - 8*k ) ) / 2 );
  size_t j = k - (2*p-1- i)*i/2;
  return(std::pair<size_t,size_t>(i,j));
}

// std::pair<size_t,size_t> row_col_t(const size_t k, const size_t p){
//   size_t i = floor( ( 2*p+1 - sqrt( (2*p+1)*(2*p+1) - 8*k ) ) / 2 );
//   size_t j = k - (2*p-1- i)*i/2;
//   return(std::pair<size_t,size_t>(i,j));
// }





tbb::concurrent_vector<typename Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex > >
sparse_LDshrink_vp(const Eigen::Map<Eigen::MatrixXd> scaled_data_a, const Eigen::Map<Eigen::MatrixXd> scaled_data_b, const Eigen::ArrayXd mapd_a, const Eigen::ArrayXd mapd_b,  const double m, const double Ne, const double cutoff,const double r2cutoff=0.01,const bool progress=false){
  using namespace tbb;
  typedef Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex> Trip;
  using trip_vec= tbb::concurrent_vector<Trip>;
  trip_vec tripletList;
  const size_t p_a = mapd_a.size();
  const size_t p_b = mapd_b.size();
  const size_t p= p_a+p_b;


  const int num_ind_a=scaled_data_a.rows();
  const int num_ind_b=scaled_data_b.rows();

  if(num_ind_a!=num_ind_b){
    Rcpp::stop("scaled_data_a and scaled_data_b must have same number of individuals (rows)");
  }
  const int num_ind=num_ind_a;
  const size_t p_od = (p_a*p_b);
  tripletList.reserve(10000);
  const bool isGeno=true;
  const double GenoMult = isGeno ? 0.5 : 1;
  const double theta=calc_theta(m);
  std::vector<double> data_vars_a(p_a);
  std::vector<double> data_vars_b(p_b);

  // Ensure that beginning of b  is greater than or equal to beginning of a
  if ((mapd_b[0]<mapd_a[0])  || (mapd_b[p_b-1]<mapd_a[p_a-1])){
    Rcpp::stop("genetic map at beginning of b must be >= map at beginning of a");
  }


  parallel_for(tbb::blocked_range<size_t>(0,p_a),[&](const blocked_range<size_t> &r){
						   for(size_t i=r.begin(); i!=r.end();i++){
						     data_vars_a[i]=(scaled_data_a.col(i).array().square().sum()/(num_ind-1))*GenoMult*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
						   }
						 });
  parallel_for(tbb::blocked_range<size_t>(0,p_b),[&](const blocked_range<size_t> &r){
						   for(size_t i=r.begin(); i!=r.end();i++){
						     data_vars_b[i]=(scaled_data_b.col(i).array().square().sum()/(num_ind-1))*GenoMult*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
						   }
						 });
  Progress pr(p_od, progress);

  parallel_for(tbb::blocked_range2d<size_t>(0,p_a,0,p_b),[&](const blocked_range2d<size_t> &r){
							   for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i ){
							     for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ) {
							       // for(size_t ir=r.begin();ir!=r.end();ir++){
							       //   auto idx =row_col_t(ir,p);
							       //   const size_t i=idx.first;
							       //   const size_t j=idx.second;
							       double map_dist=mapd_b(j)-mapd_a(i);
							       double rho = 4*Ne*(map_dist)/100;
							       rho = -rho/(2*m);
							       rho =std::exp(rho);
							       if(rho >= cutoff){
								 double cov = (scaled_data_a.col(i).dot(scaled_data_b.col(j))/(num_ind-1))*rho*GenoMult*(1-theta)*(1-theta);
								 cov=cov/std::sqrt(data_vars_a[i]*data_vars_b[j]);
								 if( ((cov*cov)>=r2cutoff) && (i!=j)){
								   tripletList.push_back(Trip(i,j,cov));
								 }
							       }
							     }
							   }
							   size_t tot_prog = (r.rows().end()-r.rows().begin())*(r.cols().end()-r.cols().begin());
							   pr.increment(tot_prog);
							 });
  return(tripletList);
}






tbb::concurrent_vector<typename Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex > >
sparse_LDshrink_v(const Eigen::Map<Eigen::MatrixXd> scaled_data, const Eigen::ArrayXd mapd,const double m, const double Ne, const double cutoff,const double r2cutoff=0.01,const bool progress=false){
  using namespace tbb;
  typedef Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex> Trip;
  using trip_vec= tbb::concurrent_vector<Trip>;
  trip_vec tripletList;
  const size_t p = mapd.size();
  const int num_ind=scaled_data.rows();
  const size_t p_od = (((p*p)-p)/2)+p;
  tripletList.reserve(10000);
  const bool isGeno=true;
  const double GenoMult = isGeno ? 0.5 : 1;
  const double theta=calc_theta(m);
  std::vector<double> data_vars(p);


  parallel_for(tbb::blocked_range<size_t>(0,p),[&](const blocked_range<size_t> &r){
  						 for(size_t i=r.begin(); i!=r.end();i++){
						   data_vars[i]=(scaled_data.col(i).array().square().sum()/(num_ind-1))*GenoMult*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
						 }
					       });
  Progress pr(p_od, progress);

  parallel_for(tbb::blocked_range<size_t>(0,p_od),[&](const blocked_range<size_t> &r){
						    for(size_t ir=r.begin();ir!=r.end();ir++){
						      auto idx =row_col_t(ir,p);
						      const size_t i=idx.first;
						      const size_t j=idx.second;
						      double map_dist=mapd(j)-mapd(i);
						      double rho = 4*Ne*(map_dist)/100;
						      rho = -rho/(2*m);
						      rho =std::exp(rho);
						      if(rho >= cutoff){
							double cov = (scaled_data.col(j).dot(scaled_data.col(i))/(num_ind-1))*rho*GenoMult*(1-theta)*(1-theta);
							cov=cov/std::sqrt(data_vars[i]*data_vars[j]);
							if( ((cov*cov)>=r2cutoff) && (i!=j)){
							  tripletList.push_back(Trip(i,j,cov));
							}
						      }
						    }
						    pr.increment(r.end()-r.begin());
						  });
  return(tripletList);
}



tbb::concurrent_vector<typename Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex > >
sparse_cor_v(const Eigen::Map<Eigen::MatrixXd> scaled_data,const double r2cutoff=0.01,const bool progress=false){
  using namespace tbb;
  typedef Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex> Trip;
  using trip_vec= tbb::concurrent_vector<Trip>;
  trip_vec tripletList;

  const size_t p = scaled_data.cols();
  const int num_ind=scaled_data.rows();
  const size_t p_od = (((p*p)-p)/2)+p;
  tripletList.reserve(10000);
  std::vector<double> data_vars(p);
  parallel_for(tbb::blocked_range<size_t>(0,p),[&](const blocked_range<size_t> &r){
  						 for(size_t i=r.begin(); i!=r.end();i++){
						   data_vars[i]=(scaled_data.col(i).array().square().sum()/(num_ind-1));
						 }
					       });
  Progress pr(p_od, progress);

  parallel_for(tbb::blocked_range<size_t>(0,p_od),[&](const blocked_range<size_t> &r){
						    for(size_t ir=r.begin();ir!=r.end();ir++){
						      auto idx =row_col_t(ir,p);
						      const size_t i=idx.first;
						      const size_t j=idx.second;
						      double cov = (scaled_data.col(j).dot(scaled_data.col(i))/(num_ind-1));
						      cov=cov/std::sqrt(data_vars[i]*data_vars[j]);
						      if( ((cov*cov)>=r2cutoff) && (i!=j)){
							tripletList.push_back(Trip(i,j,cov));
						      }
						    }
						    pr.increment(r.end()-r.begin());
						  });
  return(tripletList);
}



//[[Rcpp::export]]
SEXP sparse_LDshrink(const Eigen::Map<Eigen::MatrixXd> scaled_data, const Eigen::ArrayXd mapd,const double m, const double Ne, const double cutoff,const bool progress=false,const bool useLDshrink=true){
  const size_t p=mapd.size();
  Eigen::SparseMatrix<double> object(p,p);
  typedef Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex> Trip;
  using trip_vec= tbb::concurrent_vector<Trip>;
  trip_vec tripletList;
  if(useLDshrink){
    tripletList =  sparse_LDshrink_v(scaled_data,mapd,m,Ne,cutoff,0,progress);
  }else{
    tripletList = sparse_cor_v(scaled_data,0,progress);
  }
  for(int i=0; i<p; i++){
    tripletList.push_back({i,i,1.0});
  }
  object.setFromTriplets(tripletList.begin(), tripletList.end());
  using namespace Rcpp;
  const int    nnz = object.nonZeros();
  S4 ans("dsCMatrix");
  ans.slot("Dim")  = Dimension(object.rows(), object.cols());
            ans.slot("i") =
                IntegerVector(object.innerIndexPtr(), object.innerIndexPtr() + nnz);
            ans.slot("p")    = IntegerVector(object.outerIndexPtr(),
                                             object.outerIndexPtr() + object.outerSize() + 1);
            ans.slot("x")    = Rcpp::NumericVector(object.valuePtr(), object.valuePtr() + nnz);
  return(ans);
}



//[[Rcpp::export]]
Rcpp::DataFrame ld2df(const Eigen::Map<Eigen::MatrixXd> scaled_data, const Eigen::ArrayXd mapd,Rcpp::StringVector rsid,const double m, const double Ne, const double cutoff,const double r2cutoff=0.01,const bool progress=false,const bool useLDshrink = true){

  typedef Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex> Trip;
  using trip_vec= tbb::concurrent_vector<Trip>;
  trip_vec tripletList;
  if(useLDshrink){
    tripletList =  sparse_LDshrink_v(scaled_data,mapd,m,Ne,cutoff,r2cutoff,progress);
  }else{
    tripletList = sparse_cor_v(scaled_data,0,progress);
  }
  // auto tripletList = sparse_LDshrink_v(scaled_data,mapd,m,Ne,cutoff,r2cutoff,progress);
  const size_t num_trip=tripletList.size();
  Rcpp::StringVector rowid(num_trip);
  Rcpp::StringVector colid(num_trip);
  Rcpp::NumericVector r2vec(num_trip);


  for(size_t i=0; i<num_trip;i++){
    auto trip	= tripletList[i];
    rowid(i)=rsid(trip.row());
    colid(i)=rsid(trip.col());
    r2vec(i)=trip.value()*trip.value();
  }
  using namespace Rcpp;
  bool stringsAsFactors = false;
  return(Rcpp::DataFrame::create(_["rowsnp"]=rowid,_["colsnp"]=colid,_["r2"]=r2vec,_["stringsAsFactors"]=stringsAsFactors));
}



//[[Rcpp::export]]
Rcpp::DataFrame ld2df_p(const Eigen::Map<Eigen::MatrixXd> scaled_data_a,const Eigen::Map<Eigen::MatrixXd> scaled_data_b, const Eigen::ArrayXd mapd_a,const Eigen::ArrayXd mapd_b, Rcpp::StringVector rsid_a,Rcpp::StringVector rsid_b,const double m, const double Ne, const double cutoff,const double r2cutoff=0.01,const bool progress=false,const bool useLDshrink = true){
  typedef Eigen::Triplet<double,typename Eigen::SparseMatrix<double>::StorageIndex> Trip;
  using trip_vec= tbb::concurrent_vector<Trip>;
  trip_vec tripletList;
  if(useLDshrink){
    tripletList =  sparse_LDshrink_vp(scaled_data_a,scaled_data_b,mapd_a,mapd_b,m,Ne,cutoff,0,progress);
  }else{
    Rcpp::stop("chunkwise not implemented yet");
    tripletList = sparse_cor_v(scaled_data_a,0,progress);
  }
  // auto tripletList = sparse_LDshrink_v(scaled_data,mapd,m,Ne,cutoff,r2cutoff,progress);
  const size_t num_trip=tripletList.size();
  Rcpp::StringVector rowid(num_trip);
  Rcpp::StringVector colid(num_trip);
  Rcpp::NumericVector r2vec(num_trip);
  for(size_t i=0; i<num_trip;i++){
    auto trip	= tripletList[i];
    rowid(i)=rsid_a(trip.row());
    colid(i)=rsid_b(trip.col());
    r2vec(i)=trip.value()*trip.value();
  }
  using namespace Rcpp;
  bool stringsAsFactors = false;
  return(Rcpp::DataFrame::create(_["rowsnp"]=rowid,_["colsnp"]=colid,_["r2"]=r2vec,_["stringsAsFactors"]=stringsAsFactors));
}
