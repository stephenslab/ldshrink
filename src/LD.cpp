#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include<algorithm>
#include <functional>
#include <tuple>
#include "ldshrink.hpp"
#include <memory>
#include <progress.hpp>
//#include <mkl.h>

#include <vector>





template<typename CovFun>
LDshrinkWriter ldshrink(Genotype data,
			    const SingleList input,
			    const CovFun &cov_f,
			    const DistAnnoVec &anno,
			    const bool write_annotations=false,
			    const bool progress = false){

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

//' Internal implementation of ldshrink
//'
//'
//' @param m a number indicating the size of the panel used to create the
//genetic map ' (if using `1000-genomes-genetic-maps` from europeans, this
//number is 85)
//' @export
//[[Rcpp::export]]
SEXP ldshrink_cor(const Rcpp::NumericMatrix genotype_data,
		  const Rcpp::NumericVector anno,
                  const Rcpp::List options) {

  using namespace Rcpp;
  using index_t = int;
  const DistAnnoVec annov(anno);
  const bool isGeno = !(value_or<bool>(options,"is_haplotype", false));
  const bool write_anno = (value_or<bool>(options,"write_annotations", false));
  const bool progress = (value_or<bool>(options,"progress", false));
  std::string output_type =
    value_or<std::string>(options,"output", "data.frame");
  Genotype geno_data(genotype_data, isGeno);
  StringVector geno_names=colnames(genotype_data);
  SingleList target_snps(geno_names);
  double m = as<double>(options["m"]);
  double Ne = as<double>(options["Ne"]);
  double cutoff = as<double>(options["cutoff"]);
  LDshrinkCor cov_fun(m, Ne, cutoff, isGeno);
  LDshrinkWriter result = ldshrink(geno_data, target_snps, cov_fun,annov,write_anno,progress);
  return (result.toType(output_type));
}

//' Internal implementation of sample correlation
//'
//'
//' @param m a number indicating the size of the panel used to create the
//genetic map ' (if using `1000-genomes-genetic-maps` from europeans, this
//number is 85)
//' @export
//[[Rcpp::export]]
SEXP sample_cor(const Rcpp::NumericMatrix genotype_data,
                const Rcpp::NumericVector anno, const Rcpp::List options) {

  using namespace Rcpp;

  const DistAnnoVec annov(anno);
  const bool isGeno = !(value_or<bool>(options, "is_haplotype", false));
  const bool write_anno = (value_or<bool>(options, "write_annotations", false));
  const bool progress = (value_or<bool>(options, "progress", false));
  const std::string output_type =
      value_or<std::string>(options, "output", "data.frame");
  const double dist_cutoff =
      value_or<double>(options, "dist_cutoff", 10000000000000);
  const double cor_cutoff = value_or<double>(options, "cor_cutoff", 0);

  Genotype geno_data(genotype_data, isGeno);
  StringVector geno_names = colnames(genotype_data);
  SingleList target_snps(geno_names);

  SampleCor cov_fun(dist_cutoff, cor_cutoff);
  LDshrinkWriter result =
      ldshrink(geno_data, target_snps, cov_fun, anno, write_anno, progress);
  return (result.toType(output_type));
}
