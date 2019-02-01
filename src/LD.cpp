#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include<algorithm>
#include <functional>
#include <tuple>
#include "ldshrink.hpp"
#include <memory>
#include <progress.hpp>
#include<atomic>
#include <vector>





template<typename CovFun,typename LDWriter>
SEXP ldshrink(Genotype data,
	      const SingleList input,
	      const CovFun &cov_f,
	      const DistAnnoVec &anno,
	      const std::string output_type,
	      const bool progress = false){
  const size_t p_a= input.get_idx_row().size();
  //  auto id_pair=std::make_pair(input.get_idx_row(),input.get_idx_col());
  auto name_pair=std::make_pair(input.get_names_row(),input.get_names_col());
  const bool use_names=(name_pair.first.size()==p_a);
  LDWriter output(p_a);
  const std::size_t num_iter=input.size();
  Progress prog_bar(num_iter,progress);
  using namespace tbb;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_iter),
                    [&](const tbb::blocked_range<size_t> &r) {
		      for(size_t ir=r.begin(); ir!= r.end(); ir++){
			auto idx=input.get_idx(ir);
			const auto i=idx[0];
			const auto j=idx[1];
			if(i==j){
			  output.write_symm(idx);
			}else{
			  auto annom=anno.get(idx);
			  auto cov_m = cov_f.check(annom);
                          if (cov_m) {
                            auto res = cov_f.cor(cov_m, data.get(idx));
                            output.write_nsymm(
                                idx, std::array<double, 2>{res, annom});
                          }
                        }
                      }
                      prog_bar.increment(r.end() - r.begin());
                    });
  output.finalize();
  if (use_names) {
    return (output.toType(output_type, name_pair.first, name_pair.second));
  }else{
    return (output.toType(output_type));
  }
}

//' Internal implementation of ldshrink
//'
//'
//'
//[[Rcpp::export]]
SEXP ldshrink_cor(const Rcpp::NumericMatrix genotype_data,
		  const Rcpp::NumericVector anno,
                  const Rcpp::List options) {

  using namespace Rcpp;
  using index_t = int;
  const DistAnnoVec annov(anno);
  const bool isGeno = !(value_or<bool>(options,"is_haplotype", false));
  //  Rcpp::Rcout<<"isGeno:"<<isGeno<<std::endl;
  const bool write_anno = (value_or<bool>(options,"write_annotations", false));
  const bool progress = (value_or<bool>(options,"progress", false));
  std::string output_type =
    value_or<std::string>(options,"output", "matrix");
  Genotype geno_data(genotype_data, isGeno);
  StringVector geno_names=colnames(genotype_data);
  SingleList target_snps(geno_names);
  double m = as<double>(options["m"]);
  double Ne = as<double>(options["Ne"]);
  double cutoff = as<double>(options["cutoff"]);
  LDshrinkCor cov_fun(m, Ne, cutoff, isGeno);
  if (write_anno) {
    auto result = ldshrink<LDshrinkCor, Skyline_data_store<2>>(
        geno_data, target_snps, cov_fun, annov, output_type, progress);
    return (result);
  }
  auto result = ldshrink<LDshrinkCor, Skyline_data_store<1>>(geno_data, target_snps, cov_fun, annov, output_type, progress);
  return (result);
}

//' Internal implementation of sample correlation
//'
//'
//[[Rcpp::export]]
SEXP sample_cor(const Rcpp::NumericMatrix genotype_data,
                const Rcpp::NumericVector anno, const Rcpp::List options) {

  using namespace Rcpp;

  const DistAnnoVec annov(anno);
  const bool isGeno = !(value_or<bool>(options, "is_haplotype", false));
  const bool write_anno = (value_or<bool>(options, "write_annotations", false));
  const bool progress = (value_or<bool>(options, "progress", false));
  const std::string output_type =
      value_or<std::string>(options, "output", "matrix");
  const double dist_cutoff =
      value_or<double>(options, "dist_cutoff", 10000000000000);
  const double cor_cutoff = value_or<double>(options, "cor_cutoff", 0);

  Genotype geno_data(genotype_data, isGeno);
  StringVector geno_names = colnames(genotype_data);
  SingleList target_snps(geno_names);

  SampleCor cov_fun(dist_cutoff, cor_cutoff);
  if (write_anno) {
    auto result = ldshrink<SampleCor, Skyline_data_store<2>>(
        geno_data, target_snps, cov_fun, annov, output_type, progress);
    return (result);
  }
  auto result = ldshrink<SampleCor, Skyline_data_store<1>>(
      geno_data, target_snps, cov_fun, annov, output_type, progress);
  return (result);
}
