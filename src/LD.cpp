#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include<algorithm>
#include <functional>
#include <tuple>
#include "ldshrink.hpp"
#include <memory>
#include <progress.hpp>

#include <vector>

// [[Rcpp::depends(BH)]]


//' Calculate the constant theta, given `m`
//'
//' `calc_theta` is used to calculate the shrinkage coefficient	used by LDshrink
//'
//' @param m a number indicating the size of the panel used to create the
//genetic map ' (if using `1000-genomes-genetic-maps` from europeans, this
//number is 85)
//' @export
//[[Rcpp::export(name="calc_theta")]]
double calc_theta_exp(const double m) { return (calc_theta(m)); }

//[[Rcpp::export]]
Rcpp::NumericMatrix shrinkCov(const Rcpp::NumericMatrix S,
                              const Rcpp::NumericVector &mapd, const double m,
                              const double Ne, const double cutoff) {
  const size_t p = mapd.size();
  Rcpp::NumericMatrix nS(Rcpp::clone(S));
  Eigen::Map<Eigen::MatrixXd> mS(&nS(0, 0), p, p);
  Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> mmap(&mapd(0), p);
  ldshrinker<double> lds(mS, mmap, m, Ne, cutoff);
  lds.Shrink();
  return (nS);
}




template<typename CovFun >
LDshrinkWriter new_ldshrink(Genotype data,
			    const SingleList input,
			    const std::string anno_name,
			    const CovFun &cov_f,
			    const bool write_annotations=false,
			    const bool progress = false){




  const DistAnnoVec anno(input.get_list_elem<REALSXP>(anno_name));
  auto id_pair=std::make_pair(input.get_idx_row(),input.get_idx_col());
  auto name_pair=std::make_pair(input.get_names_row(),input.get_names_col());
  LDshrinkWriter output(id_pair,name_pair,write_annotations);
  const std::size_t num_iter=input.size();
  Progress prog_bar(num_iter,progress);
  using namespace tbb;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_iter),
                    [&](const tbb::blocked_range<size_t> &r) {
		      for(size_t ir=r.begin(); ir!= r.end(); ir++){
			auto idx=input.get_idx(ir);
			const auto i=idx.first;
			const auto j=idx.second;
			if(i==j){
			  output.write_symm(idx);
			}else{
			  auto annom=anno.get(idx);
			  auto cov_m = cov_f.check(annom);
                          if (cov_m) {
                            auto res = cov_f.cor(cov_m, data.get(idx));
                            output.write_nsymm(idx, res, annom);
                          }
                        }
                      }
		      prog_bar.increment(r.end()-r.begin());
                    });
  output.finalize();
  return (output);
}

//' Internal
//'
//'
//' @param m a number indicating the size of the panel used to create the
//genetic map ' (if using `1000-genomes-genetic-maps` from europeans, this
//number is 85)
//' @export
//[[Rcpp::export]]
SEXP ldshrink_cor(const Rcpp::List genotype_data, const Rcpp::List indices,
                  const Rcpp::List options) {

  using namespace Rcpp;
  using index_t = int;
  const bool isGeno = !(value_or<bool>(options,"is_haplotype", false));
  const bool write_anno = (value_or<bool>(options,"write_annotations", false));
  const bool progress = (value_or<bool>(options,"progress", false));
  std::string output_type =
    value_or<std::string>(options,"output", "data.frame");
  Genotype geno_data(genotype_data, indices, isGeno);
  SingleList target_snps(indices);
  double m = as<double>(options["m"]);
  double Ne = as<double>(options["Ne"]);
  double cutoff = as<double>(options["cutoff"]);
  LDshrinkCor cov_fun(m, Ne, cutoff, isGeno);
  LDshrinkWriter result = new_ldshrink(geno_data, target_snps, "map", cov_fun,write_anno,progress);
  return (result.toType(output_type));
}



//' Internal
//'
//'
//' @param m a number indicating the size of the panel used to create the
//genetic map ' (if using `1000-genomes-genetic-maps` from europeans, this
//number is 85)
//' @export
//[[Rcpp::export]]
SEXP sample_cor(const Rcpp::List genotype_data, const Rcpp::List indices,
                  const Rcpp::List options) {

  using namespace Rcpp;


    const bool isGeno = !(value_or<bool>(options,"is_haplotype", false));
  const bool write_anno = (value_or<bool>(options,"write_annotations", false));
  const std::string output_type =
    value_or<std::string>(options,"output", "data.frame");
  const double dist_cutoff   = value_or<double>(options,"dist_cutoff", 10000000000000);
  const double cor_cutoff   = value_or<double>(options,"cor_cutoff", 0);

  const bool progress = (value_or<bool>(options,"progress", false));


  Genotype geno_data(genotype_data, indices, isGeno);
  SingleList target_snps(indices);
  SampleCor cov_fun(dist_cutoff,cor_cutoff);
  LDshrinkWriter result = new_ldshrink(geno_data, target_snps, "pos", cov_fun,write_anno,progress);
  return (result.toDataFrame());
}




//[[Rcpp::export]]
Rcpp::NumericMatrix fastldshrink(const Rcpp::NumericMatrix genotype_data,
                                 const Rcpp::NumericVector &mapd,
                                 const double m, const double Ne,
                                 const double cutoff, const bool isGeno = true,
                                 const bool cov_2_cor = true) {
  const size_t p = mapd.size();
  Rcpp::NumericMatrix nS(p, p);
  Eigen::Map<Eigen::MatrixXd> mS(&nS(0, 0), p, p);
  Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> mmap(&mapd(0), p);
  ldshrinker<double> lds(mS, mmap, m, Ne, cutoff, genotype_data, isGeno);
  lds.Shrink();
  if (cov_2_cor) {
    lds.cov_2_cor();
  }
  return (nS);
}

//[[Rcpp::export]]
Eigen::MatrixXd calcDist(Eigen::VectorXd &map) {
  const size_t p = map.size();
  return (((map.transpose().colwise().replicate(p)).colwise() - map));
}
