#include "EigenH5.h"
#include <LDshrink.h>
#include<algorithm>
#include <functional>
#include <tbb/tbb.h>
#include <RcppParallel.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include<tuple>


//#include <mkl.h>
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppParallel)]]
//[[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(EigenH5)]]


typedef std::tuple<float,float ,float,float> ldp;
typedef std::tuple<int,int,int> row_col_dim;





Eigen::MatrixXd calc_cov_p(c_Matrix_external mata, c_Matrix_external matb){
  auto centereda = mata.rowwise()-mata.colwise().mean();
  auto centeredb = matb.rowwise()-matb.colwise().mean();
  return (((centereda.transpose()*centeredb)/double(mata.rows()-1)));  
}

//[[Rcpp::export(name="calc_cov_p")]]
Eigen::MatrixXd calc_cov_p_exp(Matrix_external mata, Matrix_external matb){
  auto centereda = mata.rowwise()-mata.colwise().mean();
  auto centeredb = matb.rowwise()-matb.colwise().mean();
  return (((centereda.transpose()*centeredb)/double(mata.rows()-1)));  
}


//[[Rcpp::export]]
Eigen::MatrixXd calc_cov_scaled(Matrix_external centereda, Matrix_external centeredb){
  return (((centereda.transpose()*centeredb)/double(centereda.rows()-1)));  
}





std::pair<int,int> chunk_block(const int chunk, const int chunksize, const int p){
  const int total_chunks= (chunksize+p- 1) / p;
  const int chunkstart = chunk*chunksize;
  if(chunkstart+chunksize<p){
    return(std::make_pair(chunkstart,chunksize));
  }
    return(std::make_pair(chunkstart,p-chunkstart));
}

// auto calc_chunk(const int chunksize,const int p){
//   std::vector<std::pair<int,int>> ret_vec;
//   int i=0;
//   auto nret = chunk_block(i,chunksize,p); 
//   while(nret){
//     ret_vec.push_back(nret);
//     nret=chunk_block(i++,chunksize,p);
//   }
//   return(ret_vec)
// }

constexpr int total_chunk_number(const int chunksize, const int p){
  return (chunksize+p- 1) / p;
}

auto upper_diagonal_chunks (const int chunksize, const int p){
  std::vector< std::pair< std::pair<int,int> , std::pair<int,int> > > ret_vec;
  int i=0;
  int j=0;
  const int total_chunks=total_chunk_number(chunksize,p);
  for(int i=0;i<total_chunks;i++){
    for(int j=i;j<total_chunks;j++){
      ret_vec.push_back(std::make_pair(chunk_block(i,chunksize,p),chunk_block(j,chunksize,p)));
    }
  }
  return(ret_vec);
}

void cov_2_cor_p(Matrix_internal covmat,arrayxd_internal rowvar,arrayxd_internal colvar){
//  Eigen::ArrayXd rowvar=1/covmat.diagonal().array().sqrt();
  if(covmat.cols()!=colvar.size()){
    std::cerr<<"covmat has "<<covmat.cols()<<" cols and colvar has "<<colvar.size()<<std::endl;
    std::cerr<<"covmat has "<<covmat.rows()<<" rows and rowvar has "<<rowvar.size()<<std::endl;
    Rcpp::stop("you mixed up rowvar and colvar,c:"+
      std::to_string(covmat.cols())+" "+std::to_string(colvar.size())+" r:"+
      std::to_string(covmat.rows())+" "+std::to_string(rowvar.size()));
  }
  if(covmat.rows()!=rowvar.size()){
    std::cerr<<"covmat has "<<covmat.cols()<<" cols and colvar has "<<colvar.size()<<std::endl;
    std::cerr<<"covmat has "<<covmat.rows()<<" rows and rowvar has "<<rowvar.size()<<std::endl;
    Rcpp::stop("you mixed up rowvar and colvar,c:"+
      std::to_string(covmat.cols())+" "+std::to_string(colvar.size())+" r:"+
      std::to_string(covmat.rows())+" "+std::to_string(rowvar.size()));

  }
  covmat.array().colwise()*=rowvar.sqrt().inverse();
  covmat.array().rowwise()*=colvar.sqrt().inverse().transpose();
//  covmat.diagonal().setOnes();
}

//[[Rcpp::export]]
Eigen::MatrixXd cov_2_cor_exp_p(const Eigen::MatrixXd covmat, const Eigen::ArrayXd rowvar,const Eigen::ArrayXd colvar){
  Eigen::MatrixXd retcov=covmat;
  Eigen::ArrayXd trowvar=rowvar;
  Eigen::ArrayXd tcolvar=colvar;
  cov_2_cor_p(retcov,trowvar,tcolvar);
  return(retcov);
}


//[[Rcpp::export]]
Eigen::MatrixXd eigen_dist(const Eigen::VectorXd mapa, const Eigen::VectorXd mapb){
  Eigen::MatrixXd retmat=mapa.replicate(1,mapb.size());
  return(retmat.array().rowwise()-mapb.transpose().array());
}



Eigen::MatrixXd calcLD_cov_d(const c_Matrix_external hmata,const c_arrayxd_external mapa,const c_Matrix_external hmatb,const c_arrayxd_external mapb,const arrayxd_external ldparams,const bool isDiag){
  
  
  double m,Ne,cutoff,theta;

  
  m=ldparams(0);
  Ne=ldparams(1);
  cutoff=ldparams(2);
  theta=ldparams(3);
  
  //  std::tie(Achunk,Bchunk,chunksize) = id;
  
  const int N=hmata.rows();
  const int p=mapa.size()+mapb.size();
  
  
  
  if(p!=hmata.cols()+hmatb.cols()){
    Rcpp::stop("length of mapa+mapb must be equal to number of cols of hmata+hmatb");
  }
  if(hmata.cols()!=mapa.size()){
    Rcpp::stop("length of mapa  must be equal to number of cols of hmata");
  }
  if(hmatb.cols()!=mapb.size()){
    Rcpp::stop("length of mapa  must be equal to number of cols of hmata");
  }
  
  
  double dosage_max=hmata.maxCoeff();
  bool isGeno=dosage_max>1;
  if(isGeno){
    Rcpp::stop("LD cannot currently be computed from genotype, only haplotype\n");
  }
  
  Eigen::MatrixXd S= calc_cov_p(hmata,hmatb);
  
  if(hmata.cols()!=S.rows()){
    std::cerr<<"hmata cols: "<<hmata.cols()<<" hmatb cols: "<<hmatb.cols()<<std::endl;
    std::cerr<<"S rows: "<<S.rows()<<" S cols: "<<S.cols()<<std::endl;
    Rcpp::stop("hmata.cols()!=S.rows()");
  }
  // int numSNP=S.rows();
  // if(isDiag){
  //   S.triangularView<Eigen::StrictlyLower>().setZero();
  // }
  std::vector<double> mapv;
  Eigen::ArrayXd avars(mapa.size());
  Eigen::ArrayXd bvars(mapb.size());
  
  
  
  mapv.resize(mapa.size());
  Eigen::VectorXd::Map(&mapv[0],mapa.size())=mapa;
  
  if(!is_sorted(mapv.begin(),mapv.end(),std::less<double>())){
    Rcpp::stop("Recombination map must be non-decreasing\n");
  }
  double tj=0;
  double ti=0;
  double rho=0;
  double tshrinkage;
  
int asize=mapa.size();
int bsize=mapb.size();
  if(!isDiag){
    for(int i=0; i<asize;i++){
      ti=mapa(i);
      for(int j=0; j<bsize;j++){
        tj=mapb(j);
        rho = 4*Ne*(std::abs(tj-ti))/100;
        rho=-rho/(2*m);
        tshrinkage=std::exp(rho);
        if(tshrinkage<cutoff){
          tshrinkage=0;
        }
        S(i,j)=tshrinkage*S(i,j);
      }
    }
  }else{
    for(int i=0; i<asize;i++){
      ti=mapa(i);
      for(int j=i+1; j<bsize;j++){
        tj=mapb(j);
        rho = 4*Ne*(tj-ti)/100;
        rho=-rho/(2*m);
        tshrinkage=std::exp(rho);
        if(tshrinkage<cutoff){
          tshrinkage=0;
        }
        S(i,j)=tshrinkage*S(i,j);
      }
    }
  }
  
  
  
  // // S=S+S.triangularView<Eigen::StrictlyUpper>().transpose();
  if(isDiag){
    S.triangularView<Eigen::StrictlyLower>()=S.transpose();
  }
  
  
  
  Eigen::MatrixXd SigHat = ((1-theta)*(1-theta))*S.array();
  if(isDiag){
    Eigen::ArrayXd eye(mapa.size());
    eye.setOnes();
    eye=eye.array()*(0.5*theta * (1-0.5*theta));
    SigHat.diagonal() = SigHat.diagonal().array()+eye;
//    SigHat=SigHat+eye.matrix().asDiagonal()
    // int short_dim=std::min(SigHat.rows(),SigHat.cols());
    // for(size_t i=0;i<short_dim;i++){
    //   SigHat(i,i)+=SigHat(i,i)+eye(i);
    // }
  }
  
  return(SigHat);
}





  
typedef tbb::spin_mutex Mutex;

//[[Rcpp::export]]
void calc_LD_chunk_h5(const Rcpp::DataFrame input_dff ,
                      const Rcpp::DataFrame output_dff,
                      const double m,
                      const double Ne,
                      const double cutoff,
                      const bool SNPfirst=true,
                      const bool evd=true,
                      const bool df=false,
                      const double r2cutoff=0.01){
  auto r = register_blosc(nullptr,nullptr);
  // EigenH5::start_blosc();
  using Rowmat=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
  std::unordered_map<std::string,std::shared_ptr<HighFive::File> >  m_file_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::Group> >  m_group_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::DataSet> > m_dataset_map;
  MatSlices input_f(input_dff,m_file_map,m_group_map,m_dataset_map,true);
  MatSlices output_f(output_dff,m_file_map,m_group_map,m_dataset_map,false);
  Eigen::SelfAdjointEigenSolver<Rowmat> es;
  Rowmat H;
  Rowmat S;
  Rowmat Q;
  Rowmat D;
  std::vector<std::string> rsidv;
  std::vector<std::string> rowsnp;
  std::vector<std::string> colsnp;
  std::vector<double> corv;
  std::vector<double> mapd;
  Rowmat Rsq;
  const size_t num_reg=input_f.chunk_map.size();
  Progress prog_bar(num_reg, true);

  for(auto m_it=input_f.chunk_map.begin();m_it!=input_f.chunk_map.end();m_it++){
    int chunk_id = m_it->first;
    if (Progress::check_abort() )
      Rcpp::stop("Process Aborted!");
    input_f.read_chunk(chunk_id,"dosage",H);
    if(SNPfirst){
      H.transposeInPlace();
    }
    input_f.read_chunk_vector(chunk_id,"map",mapd);
    const size_t p_chunk=H.cols();
    S.resize(p_chunk,p_chunk);
    calcLD_pa(H,mapd,S,m,Ne,cutoff);

    output_f.write_chunk(chunk_id,"R",S);
    if(df){
      input_f.read_chunk_vector(chunk_id,"SNP",rsidv);
      LD2df(S,rsidv,rowsnp,colsnp,corv,r2cutoff);
      output_f.write_chunk_vector(chunk_id,"rowSNP",rowsnp);
      output_f.write_chunk_vector(chunk_id,"colSNP",colsnp);
      output_f.write_chunk_vector(chunk_id,"r2",corv);
      rowsnp.clear();
      colsnp.clear();
      corv.clear();
    }
    if(evd){
      es.compute(S);
      Q=es.eigenvectors().rowwise().reverse();
      D=es.eigenvalues().reverse();
      Eigen::Index d_idx=0;
      Eigen::Index ta=0;
      double Dmin= D.minCoeff(&d_idx,&ta);
      if(ta!=0){
        Rcpp::stop("Indexeing has gotten screwy");
      }

      output_f.write_chunk(chunk_id,"Q",Q);
      output_f.write_chunk(chunk_id,"D",D);
    }


    Rsq=S.array().square().colwise().sum()-1;
    output_f.write_chunk(chunk_id,"L2",Rsq);
    prog_bar.increment();
  }

}
