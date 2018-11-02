#pragma once

#include <RcppEigen.h>
#include <RcppParallel.h>


template<typename T> class MatTup{
  T m_row, m_col;
  double m_value;
public:
  MatTup(const std::pair<T,T> idx):m_row(idx.first),m_col(idx.second),m_value(1.0){
  };
  MatTup(const std::pair<T,T> idx, const double value_):m_row(idx.first),m_col(idx.second),m_value(value_){
  }
  const	T row() const {return m_row;}
  const	T col() const {return m_col;}
  const	double value() const {return m_value;}
  const std::pair<T,T> idx() const {return std::make_pair(m_row, m_col);}
};

// template<int N> class AnnoTup{
//   using T=int;
//   T m_row, m_col;
//   std::array<double,N> m_value;
// public:
//   AnnoTup(const std::pair<T, T> &idx)
//     : m_row(idx.first), m_col(idx.second){
//     m_value.fill(1.0);
//   };
//   AnnoTup(const std::pair<T, T> &idx,const std::array<double,N> data)
//     : m_row(idx.first), m_col(idx.second), m_value(data) {}
//   const	T row() const {return m_row;}
//   const	T col() const {return m_col;}
//   const double value() const { return m_value[0];}
//   const	std::array<double,N-1> anno_value() const {
//     std::array<double,N-1> ret;
//     std::copy(m_value.begin()+1,m_value.end(),ret.begin());
//     return ret;
//   }

// };



class LDshrinkWriter{
protected:
  using T=int;
  tbb::concurrent_vector<MatTup<T> > data_store;
  tbb::concurrent_vector<MatTup<T> > anno_data_store;
  using	name_vec=Rcpp::StringVector;
  const name_vec snpnames_a;
  const name_vec snpnames_b;
  const size_t p_a;
  const size_t p_b;
  const bool write_annotations;
public:
  LDshrinkWriter(const std::vector<T> &inp,const name_vec names,const bool write_annotations=false);
  LDshrinkWriter(const std::pair<std::vector<T>,std::vector<T>> &inpconst, const std::pair<name_vec,name_vec> names,const bool write_annotations=false);
  void write_symm(const std::pair<T,T> idx);
  void write_nsymm(const std::pair<T,T> idx,const double res,const double annom);
  Eigen::SparseMatrix<double> sparseMatrix() const;
  void finalize();
  Rcpp::S4 dsCMatrix() const;
  Rcpp::DataFrame toDataFrame(const bool use_rownames=true)const;
  Rcpp::DataFrame toAnnotationDataFrame(const bool use_rownames)const;
  SEXP toType(const std::string	output_type)const;
};

