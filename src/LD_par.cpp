#include "EigenH5.h"
#include <LDshrink.h>
#include<algorithm>
#include <functional>
#include <progress.hpp>
#include <progress_bar.hpp>
#include<tuple>
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(EigenH5)]]




//[[Rcpp::export]]
void setLDshrinkThreads(int n){
#ifdef _OPENMP
  Eigen::setNbThreads(n);
#else
  Rcpp::Rcerr<<"LDshrink compiled without openmp support, only one thread will be used"<<std::endl;
#endif
}

//[[Rcpp::export]]
int getLDshrinkThreads(){
  return(Eigen::nbThreads( ));
}









// void calc_LD_chunk_h5(const Rcpp::DataFrame input_dff ,
//                       const Rcpp::DataFrame output_dff,
//                       const double m,
//                       const double Ne,
//                       const double cutoff,
//                       const bool SNPfirst=true,
//                       const bool evd=true,
//                       const bool svd=false,
//                       const bool df=false,
//                       const double r2cutoff=0.01){
//   auto r = register_blosc(nullptr,nullptr);

//   IO_h5 IO_obj(input_dff,output_dff);

  
//   int nt= getLDshrinkThreads();
//   Rcpp::Rcerr<<"Using :"<<nt<<" threads"<<std::endl;

//   const size_t num_reg=IO_obj.input_f.chunk_map.size();
//   int mm;

//   using mmat=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
//   using marray=Eigen::ArrayXd;
//   using mmati=Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>;


//   Progress prog_bar(num_reg, true);
//   Eigen::MatrixXd H;
//   Eigen::MatrixXd Rsq;
//   Eigen::MatrixXd S;
//   Eigen::MatrixXd Q;
//   Eigen::MatrixXd D;
//   std::vector<double> mapd;
//   std::vector<std::string> rsidv;
//   std::vector<std::string> rowsnp;
//   std::vector<std::string> colsnp;
//   std::vector<double> corv;
//   Eigen::SelfAdjointEigenSolver<mmat> es;
//   Eigen::JacobiSVD<mmat> svdc;
//   LDshrinker<double> lds(m,Ne,cutoff);


//   for(auto m_it=IO_obj.input_f.chunk_map.begin();m_it!=IO_obj.input_f.chunk_map.end();m_it++){
//     int chunk_id = m_it->first;
//     if (Progress::check_abort() )
//       Rcpp::stop("Process Aborted!");
//     IO_obj.input_f.read_chunk(chunk_id,"dosage",H);
//     if(SNPfirst){
//       H.transposeInPlace();
//     }
//     IO_obj.input_f.read_chunk_vector(chunk_id,"map",mapd);
//     const size_t p_chunk=H.cols();
//     const int n=H.rows();
//     S.resize(p_chunk,p_chunk);
//     S=LDshrinker<double>::calc_cov(H);




//     if(svd){
//       mmat tH=H.eval();
//       svdc.compute(tH,Eigen::ComputeThinV);
//       D=svdc.singularValues();
//       Q=svdc.matrixV();
//       Q.transposeInPlace();
//       IO_obj.output_f.write_chunk(chunk_id,"V",Q);
//       IO_obj.output_f.write_chunk(chunk_id,"d",D);
//     }

//     IO_obj.output_f.write_chunk(chunk_id,"R",S);
//     Rsq=S.array().square().colwise().sum()-1;
//     IO_obj.output_f.write_chunk(chunk_id,"L2",Rsq);
    
//     if(df){
//       IO_obj.input_f.read_chunk_vector(chunk_id,"SNP",rsidv);
//       LD2df(S,rsidv,rowsnp,colsnp,corv,r2cutoff);
//       IO_obj.output_f.write_chunk_vector(chunk_id,"rowSNP",rowsnp);
//       IO_obj.output_f.write_chunk_vector(chunk_id,"colSNP",colsnp);
//       IO_obj.output_f.write_chunk_vector(chunk_id,"r2",corv);
//       rowsnp.clear();
//       colsnp.clear();
//       corv.clear();
//     }
//     if(evd){
      
//       es.compute(S);
//       auto res = es.info();
//       Q=es.eigenvectors();
//       D=es.eigenvalues();
//       Q.rowwise().reverseInPlace();
//       D.reverseInPlace();
//       Eigen::Index d_idx=0;
//       Eigen::Index ta=0;
      
      
//       double Dmin= D.minCoeff(&d_idx,&ta);
//       if(ta!=0){
//         Rcpp::stop("Indexeing has gotten screwy");
//       }

//       if(res!=Eigen::Success){
//         Rcpp::stop("EVD unsuccessful!");
//       }

//       IO_obj.output_f.write_chunk(chunk_id,"Q",Q);
//       IO_obj.output_f.write_chunk(chunk_id,"D",D);
//     }
//     prog_bar.increment();
//   }

// }
