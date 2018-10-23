// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_theta_exp
double calc_theta_exp(const double m);
RcppExport SEXP _ldshrink_calc_theta_exp(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_theta_exp(m));
    return rcpp_result_gen;
END_RCPP
}
// shrinkCov
Rcpp::NumericMatrix shrinkCov(const Rcpp::NumericMatrix S, const Rcpp::NumericVector& mapd, const double m, const double Ne, const double cutoff);
RcppExport SEXP _ldshrink_shrinkCov(SEXP SSEXP, SEXP mapdSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mapd(mapdSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(shrinkCov(S, mapd, m, Ne, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// fastldshrink
Rcpp::NumericMatrix fastldshrink(const Rcpp::NumericMatrix genotype_data, const Rcpp::NumericVector& mapd, const double m, const double Ne, const double cutoff, const bool isGeno, const bool cov_2_cor);
RcppExport SEXP _ldshrink_fastldshrink(SEXP genotype_dataSEXP, SEXP mapdSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP, SEXP isGenoSEXP, SEXP cov_2_corSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type genotype_data(genotype_dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mapd(mapdSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const bool >::type isGeno(isGenoSEXP);
    Rcpp::traits::input_parameter< const bool >::type cov_2_cor(cov_2_corSEXP);
    rcpp_result_gen = Rcpp::wrap(fastldshrink(genotype_data, mapd, m, Ne, cutoff, isGeno, cov_2_cor));
    return rcpp_result_gen;
END_RCPP
}
// calcDist
Eigen::MatrixXd calcDist(Eigen::ArrayXd& map);
RcppExport SEXP _ldshrink_calcDist(SEXP mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type map(mapSEXP);
    rcpp_result_gen = Rcpp::wrap(calcDist(map));
    return rcpp_result_gen;
END_RCPP
}
// sparse_ldshrink
SEXP sparse_ldshrink(Eigen::MatrixXd data, std::vector<double> mapd, Rcpp::IntegerVector indices, const double m, const double Ne, const double cutoff, const int total_size, const bool progress, const bool useldshrink);
RcppExport SEXP _ldshrink_sparse_ldshrink(SEXP dataSEXP, SEXP mapdSEXP, SEXP indicesSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP, SEXP total_sizeSEXP, SEXP progressSEXP, SEXP useldshrinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type mapd(mapdSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const int >::type total_size(total_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< const bool >::type useldshrink(useldshrinkSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_ldshrink(data, mapd, indices, m, Ne, cutoff, total_size, progress, useldshrink));
    return rcpp_result_gen;
END_RCPP
}
// sparse_ldshrink_p
SEXP sparse_ldshrink_p(Eigen::MatrixXd data_a, Eigen::MatrixXd data_b, std::vector<double> mapd_a, std::vector<double> mapd_b, Rcpp::IntegerVector indices_a, Rcpp::IntegerVector indices_b, const double m, const double Ne, const double cutoff, const int total_size, const bool progress, const bool useldshrink);
RcppExport SEXP _ldshrink_sparse_ldshrink_p(SEXP data_aSEXP, SEXP data_bSEXP, SEXP mapd_aSEXP, SEXP mapd_bSEXP, SEXP indices_aSEXP, SEXP indices_bSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP, SEXP total_sizeSEXP, SEXP progressSEXP, SEXP useldshrinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data_a(data_aSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data_b(data_bSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type mapd_a(mapd_aSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type mapd_b(mapd_bSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type indices_a(indices_aSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type indices_b(indices_bSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const int >::type total_size(total_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< const bool >::type useldshrink(useldshrinkSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_ldshrink_p(data_a, data_b, mapd_a, mapd_b, indices_a, indices_b, m, Ne, cutoff, total_size, progress, useldshrink));
    return rcpp_result_gen;
END_RCPP
}
// ld2df
Rcpp::DataFrame ld2df(Eigen::MatrixXd data, std::vector<double> mapd, Rcpp::RObject rsid, const double m, const double Ne, const double cutoff, const double r2cutoff, const bool progress, const bool useldshrink);
RcppExport SEXP _ldshrink_ld2df(SEXP dataSEXP, SEXP mapdSEXP, SEXP rsidSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP, SEXP r2cutoffSEXP, SEXP progressSEXP, SEXP useldshrinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type mapd(mapdSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type rsid(rsidSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const double >::type r2cutoff(r2cutoffSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< const bool >::type useldshrink(useldshrinkSEXP);
    rcpp_result_gen = Rcpp::wrap(ld2df(data, mapd, rsid, m, Ne, cutoff, r2cutoff, progress, useldshrink));
    return rcpp_result_gen;
END_RCPP
}
// ld2df_p
Rcpp::DataFrame ld2df_p(Eigen::MatrixXd data_a, Eigen::MatrixXd data_b, std::vector<double> mapd_a, std::vector<double> mapd_b, Rcpp::RObject rsid_a, Rcpp::RObject rsid_b, const double m, const double Ne, const double cutoff, const double r2cutoff, const bool progress, const bool useldshrink);
RcppExport SEXP _ldshrink_ld2df_p(SEXP data_aSEXP, SEXP data_bSEXP, SEXP mapd_aSEXP, SEXP mapd_bSEXP, SEXP rsid_aSEXP, SEXP rsid_bSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP, SEXP r2cutoffSEXP, SEXP progressSEXP, SEXP useldshrinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data_a(data_aSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data_b(data_bSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type mapd_a(mapd_aSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type mapd_b(mapd_bSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type rsid_a(rsid_aSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type rsid_b(rsid_bSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const double >::type r2cutoff(r2cutoffSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< const bool >::type useldshrink(useldshrinkSEXP);
    rcpp_result_gen = Rcpp::wrap(ld2df_p(data_a, data_b, mapd_a, mapd_b, rsid_a, rsid_b, m, Ne, cutoff, r2cutoff, progress, useldshrink));
    return rcpp_result_gen;
END_RCPP
}
// flip_allele
Rcpp::LogicalVector flip_allele(const Rcpp::IntegerVector& gwas_ref, const Rcpp::IntegerVector& gwas_alt, const Rcpp::IntegerVector& ld_ref, const Rcpp::IntegerVector& ld_alt);
RcppExport SEXP _ldshrink_flip_allele(SEXP gwas_refSEXP, SEXP gwas_altSEXP, SEXP ld_refSEXP, SEXP ld_altSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type gwas_ref(gwas_refSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type gwas_alt(gwas_altSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ld_ref(ld_refSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ld_alt(ld_altSEXP);
    rcpp_result_gen = Rcpp::wrap(flip_allele(gwas_ref, gwas_alt, ld_ref, ld_alt));
    return rcpp_result_gen;
END_RCPP
}
// sorted_snp_df
bool sorted_snp_df(const Rcpp::DataFrame& snp_info);
RcppExport SEXP _ldshrink_sorted_snp_df(SEXP snp_infoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type snp_info(snp_infoSEXP);
    rcpp_result_gen = Rcpp::wrap(sorted_snp_df(snp_info));
    return rcpp_result_gen;
END_RCPP
}
// set_ld_region
Rcpp::IntegerVector set_ld_region(const Rcpp::DataFrame& ld_regions, const Rcpp::DataFrame& snp_info, const bool assign_all);
RcppExport SEXP _ldshrink_set_ld_region(SEXP ld_regionsSEXP, SEXP snp_infoSEXP, SEXP assign_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type ld_regions(ld_regionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type snp_info(snp_infoSEXP);
    Rcpp::traits::input_parameter< const bool >::type assign_all(assign_allSEXP);
    rcpp_result_gen = Rcpp::wrap(set_ld_region(ld_regions, snp_info, assign_all));
    return rcpp_result_gen;
END_RCPP
}
// interpolate_map
Rcpp::NumericVector interpolate_map(const Rcpp::NumericVector& map, const Rcpp::IntegerVector map_pos, const Rcpp::IntegerVector target_pos, const bool progress);
RcppExport SEXP _ldshrink_interpolate_map(SEXP mapSEXP, SEXP map_posSEXP, SEXP target_posSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type map_pos(map_posSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type target_pos(target_posSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(interpolate_map(map, map_pos, target_pos, progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ldshrink_calc_theta_exp", (DL_FUNC) &_ldshrink_calc_theta_exp, 1},
    {"_ldshrink_shrinkCov", (DL_FUNC) &_ldshrink_shrinkCov, 5},
    {"_ldshrink_fastldshrink", (DL_FUNC) &_ldshrink_fastldshrink, 7},
    {"_ldshrink_calcDist", (DL_FUNC) &_ldshrink_calcDist, 1},
    {"_ldshrink_sparse_ldshrink", (DL_FUNC) &_ldshrink_sparse_ldshrink, 9},
    {"_ldshrink_sparse_ldshrink_p", (DL_FUNC) &_ldshrink_sparse_ldshrink_p, 12},
    {"_ldshrink_ld2df", (DL_FUNC) &_ldshrink_ld2df, 9},
    {"_ldshrink_ld2df_p", (DL_FUNC) &_ldshrink_ld2df_p, 12},
    {"_ldshrink_flip_allele", (DL_FUNC) &_ldshrink_flip_allele, 4},
    {"_ldshrink_sorted_snp_df", (DL_FUNC) &_ldshrink_sorted_snp_df, 1},
    {"_ldshrink_set_ld_region", (DL_FUNC) &_ldshrink_set_ld_region, 3},
    {"_ldshrink_interpolate_map", (DL_FUNC) &_ldshrink_interpolate_map, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ldshrink(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
