#include <range/v3/all.hpp> 
#include <LDshrink.h>
#include<algorithm>
#include <functional>
#include <tbb/tbb.h>
#include <RcppParallel.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include<tuple>
#include <EigenH5.h>


//#include <mkl.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]

typedef std::tuple<float,float ,float,float> ldp;
typedef std::tuple<int,int,int> row_col_dim;

//[[Rcpp::export]]
Rcpp::IntegerVector test_cumsum(Rcpp::IntegerVector ld_region){
  using namespace ranges;
  
  // iterator_range<Rcpp::IntegerVector::iterator> rgi {ld_region.begin(),ld_region.end()};

  
  iterator_range<Rcpp::IntegerVector::iterator> rgi {ld_region.begin(),ld_region.end()};
  int cnt = ld_region[0];
  std::vector<int> retvec = view::partial_sum(rgi, [&cnt](int i, int j) { 
    Rcpp::Rcout<<"cnt was: "<<cnt<<std::endl;
    cnt+=j;
    Rcpp::Rcout<<"cnt is: "<<cnt<<std::endl;
    
    return(cnt);
  });
  
  
  return(Rcpp::wrap(retvec));
}


//[[Rcpp::export]]
Rcpp::IntegerMatrix ld_chunks(Rcpp::IntegerVector ld_region){

  using namespace ranges;
  iterator_range<Rcpp::IntegerVector::iterator> ld_r {ld_region.begin(),ld_region.end()};
  //std::vector<int> ld_r(ld_region.begin(),ld_region.end());
  auto rng= view::zip(view::ints(0),ld_r) | view::group_by([](auto a, auto b) {
    return std::get<1>(a) == std::get<1>(b);});
  auto retr=view::transform(rng,[](auto r){
      int g0 = std::get<0>(ranges::index(r,0));
      return(view::concat(view::single(g0),view::single(distance(r))));});
  
  std::vector<int> retv = view::join(retr);
  const size_t n_ld=retv.size()/2;
  Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >tm(retv.data(),n_ld,2);
  Rcpp::IntegerMatrix retm= Rcpp::wrap(tm);
 
  return(Rcpp::wrap(retm));  
  //  std::vector<int>
}


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

typedef tbb::spin_mutex Mutex;



//[[Rcpp::export]]
void calc_ld_h5_exp(const std::string input_file,
                    const std::string output_file,
		    const std::vector<int> &ld_region,
		    const std::vector<double> &mapd,
                    const double m,
		    const double Ne,
		    const double cutoff,
		    const bool SNPfirst=true){

  using namespace HighFive;
  
  auto Xf=File(input_file,File::ReadOnly);
  auto Sf=File(output_file,File::ReadWrite|File::Create);
  
  auto Xd=Xf.getDataSet("dosage");
  auto X_dim = Xd.getDataDimensions();
  const size_t N= SNPfirst ? X_dim[1] : X_dim[0];
  const size_t p=SNPfirst ? X_dim[0] : X_dim[1];
  
  //  auto Md=Xf.getGroup("SNPinfo").getDataSet("map");
  auto r = register_blosc(nullptr,nullptr);
  //  std::vector<int> ld_region;
  //  Xf.getGroup("SNPinfo").getDataSet("region_id").read(ld_region);
  if(ld_region.size()!=p){
    Rcpp::Rcerr<<"Size of ld_region: "<<ld_region.size()<<std::endl;
    Rcpp::Rcerr<<"SNPs in dosage: "<<p<<std::endl;
    Rcpp::stop("size of SNPinfo/region_ids is not equal to number of SNPs in dosage");
  }
  
  
  using namespace ranges;
  //iterator_range<Rcpp::IntegerVector::iterator> ld_r {ld_region.begin(),ld_region.end()};
  //std::vector<int> ld_r(ld_region.begin(),ld_region.end());
  Mutex mutex;
  auto ldshrink_lambda_evd = [&](const std::string& sgname,const size_t begin, const size_t chunk_size){
    std::vector<double> tmap(chunk_size);
    std::copy_n(mapd.begin()+begin,chunk_size,tmap.begin());
    Eigen::MatrixXd X;
    Eigen::MatrixXd S(chunk_size,chunk_size);
    {
      tbb::spin_mutex::scoped_lock lock(mutex);
      if(SNPfirst){
        Xd.selectEigen({begin,0},{chunk_size,N},{}).read(X);
        X.transposeInPlace();
      }else{
        Xd.selectEigen({0,begin},{N,chunk_size},{}).read(X);
      }
      //      Md.select({begin},{chunk_size},{}).read(map);
    }
    
    calcLD_pa(X,tmap,S,m,Ne,cutoff);
    Eigen::VectorXd Rsq=S.array().square().colwise().sum()-1;
    {      
      tbb::spin_mutex::scoped_lock lock(mutex);
      Filter filter({chunk_size,chunk_size},FILTER_BLOSC,1);
      
      auto Sd=Sf.createOrGetGroup("LD").createOrGetGroup(sgname).createDataSet("R",DataSpace::From(S),AtomicType<double>(),filter.getId(),false);
      Sd.write(S);
      
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    //Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(S);
    Eigen::MatrixXd Q=es.eigenvectors().rowwise().reverse();
    Eigen::VectorXd D=es.eigenvalues().reverse();
    Eigen::Index d_idx=0;
    double Dmin= D.minCoeff(&d_idx);
    if(Dmin<=0){
      Rcpp::Rcerr<<"In ld_region "<<sgname<<"from "<<begin<<"to "<<begin+chunk_size+1<<std::endl;
      Rcpp::Rcerr<<"Encountered non-positive eigenvalue: "<<Dmin<<" at index: "<<d_idx<<std::endl;
      Rcpp::stop("Min of eigenvalue must be positive");
    }
    {      
      tbb::spin_mutex::scoped_lock lock(mutex);
      Filter filter({chunk_size,chunk_size},FILTER_BLOSC,1);
      auto Qd=Sf.createOrGetGroup("EVD").createOrGetGroup(sgname).createDataSet("Q",DataSpace::From(S),AtomicType<double>(),filter.getId(),false);
      Qd.write(Q);
      Filter filterd({chunk_size},FILTER_BLOSC,1);
      auto Dd=Sf.createOrGetGroup("EVD").createOrGetGroup(sgname).createDataSet("D",DataSpace::From(D),AtomicType<double>(),filterd.getId(),false);
      Dd.write(D);
      auto Rsqd=Sf.createOrGetGroup("LDSC").createOrGetGroup(sgname).createDataSet("L2",DataSpace::From(Rsq),AtomicType<double>(),filterd.getId(),false);
      Rsqd.write(Rsq);
    }
    
    
  };
  
  auto rng= view::zip(view::ints(0),ld_region) | view::group_by([](auto a, auto b) {
      return std::get<1>(a) == std::get<1>(b);});
  
  const size_t num_reg = distance(rng);
  Progress prog_bar(num_reg, true);
  ranges::for_each(rng,[&](auto r){
    auto first_el=ranges::index(r,0);
    const size_t offset = static_cast<size_t>(std::get<0>(first_el));
    const int tld_region = std::get<1>(first_el);
    if(!Rcpp::IntegerVector::is_na(tld_region)){
      const size_t chunksize = distance(r);
      ldshrink_lambda_evd(std::to_string(tld_region),offset,chunksize);
      prog_bar.increment();
    }
    });
}

// 
// class LDpar {
//   const int N;
//   const int p;
//   const double theta;
//   const arrayxd_external ldparams;
//   const double chunksize;
//   const std::string out_file;
//   const std::string group_name;
//   const std::string data_name;
//   mutable Mutex* mutex;
//   const std::vector< std::pair< std::pair<int,int> , std::pair<int,int> > > chunks_vec;
//   public:
//   LDpar(const arrayxd_external _ldparams,const std::vector<std::pair<size_t,size_t> > ld_chunks,const std::vector<std::string> names,Mutex &ref_mut):
//     N(hmat.rows()),
//     p(hmat.cols()),
//     ldparams(_ldparams),
//     theta(_ldparams(3)),
//     chunksize(_chunksize),
//     out_file(names[0]),
//     group_name(names[1]),
//     data_name(names[2]),
//     mutex(&ref_mut),
//     chunks_vec(upper_diagonal_chunks(chunksize,p))
//     {
//     
//   }
//     void operator()(const tbb::blocked_range<size_t> &r) const{
//       for(size_t i=r.begin();i!=r.end();i++){
//         auto p_pairs=chunks_vec[i];
//         unsigned long chunkstart_row = p_pairs.first.first;
//         unsigned long chunksize_row = p_pairs.first.second;
//         unsigned long chunkstart_col = p_pairs.second.first;
//         unsigned long chunksize_col = p_pairs.second.second;
//         const auto Ablock=chunk_block(chunkstart_row,chunksize,p);
//         const auto Bblock=chunk_block(chunkstart_col,chunksize,p);
//         const bool isDiag=chunkstart_row==chunkstart_col;
//         
//         
//         const c_Matrix_external hmata(hmat.block(0,Ablock.first,N,Ablock.second).data(),N,Ablock.second);
//         const c_Matrix_external hmatb(hmat.block(0,Bblock.first,N,Bblock.second).data(),N,Bblock.second);
//         
//         const c_arrayxd_external mapa(map.segment(Ablock.first,Ablock.second).data(),Ablock.second);
//         const c_arrayxd_external mapb(map.segment(Bblock.first,Bblock.second).data(),Bblock.second);
//         Eigen::MatrixXd covmat =calcLD_cov_d(hmata, mapa,  hmatb, mapb, ldparams, isDiag);
//         if(isDiag){
//           Eigen::ArrayXd covdiag= covmat.diagonal();
//           cov_2_cor_p(covmat,covdiag,covdiag);
//           covmat.diagonal().setOnes();
//         }else{
//           Eigen::ArrayXd cov_a=calc_variance(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
//           Eigen::ArrayXd cov_b=calc_variance(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
//           cov_2_cor_p(covmat,cov_a,cov_b);
//         }
//         {
//           using namespace HighFive;
//           tbb::spin_mutex::scoped_lock lock(*mutex);
//           File file(out_file,File::ReadWrite);
//           auto grp = file.getGroup(group_name);
//           auto dataset = grp.getDataSet(data_name);
//           dataset.select({chunkstart_row,chunkstart_col},{chunksize_row,chunksize_col}).write(covmat);
//           if(!isDiag){
//             // covmat.transposeInPlace();
//             // dataset.select({chunkstart_col,chunkstart_row},{chunksize_col,chunksize_row}).write(covmat);
//           }
//       }
//     }
//     }
//     
// };



// 
// 
// 
// void calcLD_par_h5(const Matrix_external hmat,const arrayxd_external map,const arrayxd_external ldparams,const int chunksize,const std::vector<std::string> grp_names){
//   
//   auto  out_file = grp_names[0];
//   auto group_name = grp_names[1];
//   auto data_name = grp_names[2];
//   
//   const int N=hmat.rows();
//   const int p=map.size();
//   
//   using namespace HighFive;
//   {
//   File file(out_file,File::ReadWrite,File::Create);
//   auto grp = file.createGroup(group_name);
//   
//   // char *version, *date;
//   // int r;
//   
//   /* Register the filter with the library */
//   // r = register_blosc(&version, &date);
//   // std::vector<size_t> cshape{1000,1000};
//   hsize_t chunk[2] = {1000,1000};
//   auto dcpl = H5Pcreate(H5P_DATASET_CREATE);
//   auto status = H5Pset_deflate (dcpl, 9);
//   status = H5Pset_chunk (dcpl, 2, chunk);
//   // Filter filter(cshape,H5Z_FILTER_DEFLATE,r);
//   auto dataset = grp.createDataSet(data_name,DataSpace({(unsigned long)p,(unsigned long)p}),AtomicType<double>(),dcpl);
//   }
//   auto chunks_vec=upper_diagonal_chunks(chunksize,p);
//   size_t num_chunks_tot = chunks_vec.size();
//   tbb::spin_mutex mutex;
//   LDpar ldp(hmat,map,ldparams,chunksize,grp_names,mutex);
//   tbb::parallel_for(tbb::blocked_range<size_t>(0,num_chunks_tot), ldp);
// 
// }

// 
// 
// Eigen::MatrixXd calcLD_par(const Matrix_external hmat,const arrayxd_external map,const arrayxd_external ldparams,const arrayxi_external id){
// 
// 
//   // Need to decide whether it's worth it to scale the matrices ahead of time or not (if we do, it's a little harder to tell if they're genotype or haplotype)
//   double m,Ne,cutoff,theta;
//   int Achunk,Bchunk,chunksize;
//   
//   m=ldparams(0);
//   Ne=ldparams(1);
//   cutoff=ldparams(2);
//   theta=ldparams(3);
//   
//   Achunk=id(0);
//   Bchunk=id(1);
//   chunksize=id(2);
//   //  std::tie(Achunk,Bchunk,chunksize) = id;
//   
//   const int N=hmat.rows();
//   const int p=map.size();
//   
//   
//   
//   if(p!=hmat.cols()){
//     Rcpp::stop("length of map must be equal to number of cols of hmat");
//   }
//   
//   
//   const auto Ablock=chunk_block(Achunk,chunksize,p);
//   const auto Bblock=chunk_block(Bchunk,chunksize,p);
//   const bool isDiag=Achunk==Bchunk;
// 
// 
//   const c_Matrix_external hmata(hmat.block(0,Ablock.first,N,Ablock.second).data(),N,Ablock.second);
//   const c_Matrix_external hmatb(hmat.block(0,Bblock.first,N,Bblock.second).data(),N,Bblock.second);
//   
//   const c_arrayxd_external mapa(map.segment(Ablock.first,Ablock.second).data(),Ablock.second);
//   const c_arrayxd_external mapb(map.segment(Bblock.first,Bblock.second).data(),Bblock.second);
//   Eigen::MatrixXd covmat =calcLD_cov_d(hmata, mapa,  hmatb, mapb, ldparams, isDiag);
//   if(isDiag){
//     Eigen::ArrayXd covdiag= covmat.diagonal();
//     cov_2_cor_p(covmat,covdiag,covdiag);
//     covmat.diagonal().setOnes();
//   }else{
//     Eigen::ArrayXd cov_a=calc_variance(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
//     Eigen::ArrayXd cov_b=calc_variance(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
//     cov_2_cor_p(covmat,cov_a,cov_b);
//   }
//   return(covmat);
// }
