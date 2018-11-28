#pragma once
#include <Rcpp.h>

template<typename D>
class LDshrinkWriter{
protected:
  using T=int;
  D data_store;
  using	name_vec=Rcpp::StringVector;
  const name_vec snpnames_row;
  const name_vec snpnames_col;
  const size_t p_row;
  const size_t p_col;
  const bool write_annotations;
public:
  LDshrinkWriter(const std::vector<T> &inp,const name_vec names,const bool write_annotations=false);
  LDshrinkWriter(const std::pair<std::vector<T>,std::vector<T>> &inpconst, const std::pair<name_vec,name_vec> names,const bool write_annotations=false);
  void write_symm(const std::pair<T,T> idx);
  void write_nsymm(const std::pair<T,T> idx,const double res,const double annom);
  Eigen::SparseMatrix<double> sparseMatrix() const;
  void finalize();
  SEXP dsCMatrix() const;
  Rcpp::DataFrame toDataFrame(const bool use_rownames=true)const;
  Rcpp::DataFrame toAnnotationDataFrame(const bool use_rownames)const;
  SEXP toType(const std::string	output_type)const;
};
