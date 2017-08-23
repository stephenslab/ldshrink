#include <LDshrink.h>
#include<algorithm>
#include <functional>
#include <mkl.h>
// [[Rcpp::depends(RcppEigen)]]

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

Eigen::ArrayXd calc_variance(c_Matrix_external mat){
  int n=mat.rows();
  return(((mat.rowwise()-(mat.colwise().mean())).array().square().colwise().sum())/(n-1));
}


std::pair<int,int> chunk_block(const int chunk, const int chunksize, const int p){
  const int chunkstart = chunk*chunksize;
  if(chunkstart+chunksize<p){
    return(std::make_pair(chunkstart,chunksize));
  }else{
    return(std::make_pair(chunkstart,p-chunkstart));
  }
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
  
  
  
  // avars=calc_variance(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  // bvars=calc_variance(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  // if(avars.size()!=mapa.size()){
  //   Rcpp::stop("avars.size != mapa.size()");
  // }
  // if(bvars.size()!=mapb.size()){
  //   Rcpp::stop("bvars.size != mapb.size()");
  // }
  
  // //       % SigHat is derived from Li and Stephens model (2003)
  // //         SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  // // SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
  // cov_2_cor_p(SigHat,avars,bvars);
  // // if(isDiag){
  // //   SigHat.diagonal().setOnes();
  // // }
  return(SigHat);
}



//[[Rcpp::export]]
Eigen::MatrixXd calcLD_par(const Matrix_external hmat,const arrayxd_external map,const arrayxd_external ldparams,const arrayxi_external id){


  // Need to decide whether it's worth it to scale the matrices ahead of time or not (if we do, it's a little harder to tell if they're genotype or haplotype)
  double m,Ne,cutoff,theta;
  int Achunk,Bchunk,chunksize;
  
  m=ldparams(0);
  Ne=ldparams(1);
  cutoff=ldparams(2);
  theta=ldparams(3);
  
  Achunk=id(0);
  Bchunk=id(1);
  chunksize=id(2);
  //  std::tie(Achunk,Bchunk,chunksize) = id;
  
  const int N=hmat.rows();
  const int p=map.size();
  
  
  
  if(p!=hmat.cols()){
    Rcpp::stop("length of map must be equal to number of cols of hmat");
  }
  
  
  const auto Ablock=chunk_block(Achunk,chunksize,p);
  const auto Bblock=chunk_block(Bchunk,chunksize,p);
  const bool isDiag=Achunk==Bchunk;


  const c_Matrix_external hmata(hmat.block(0,Ablock.first,N,Ablock.second).data(),N,Ablock.second);
  const c_Matrix_external hmatb(hmat.block(0,Bblock.first,N,Bblock.second).data(),N,Bblock.second);
  
  const c_arrayxd_external mapa(map.segment(Ablock.first,Ablock.second).data(),Ablock.second);
  const c_arrayxd_external mapb(map.segment(Bblock.first,Bblock.second).data(),Bblock.second);
  Eigen::MatrixXd covmat =calcLD_cov_d(hmata, mapa,  hmatb, mapb, ldparams, isDiag);
  if(isDiag){
    Eigen::ArrayXd covdiag= covmat.diagonal();
    cov_2_cor_p(covmat,covdiag,covdiag);
    covmat.diagonal().setOnes();
  }else{
    Eigen::ArrayXd cov_a=calc_variance(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
    Eigen::ArrayXd cov_b=calc_variance(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
    cov_2_cor_p(covmat,cov_a,cov_b);
  }
  return(covmat);
}
