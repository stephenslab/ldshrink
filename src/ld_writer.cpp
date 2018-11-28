// #pragma once
// #include <ldshrink.hpp>
// #include "tbb/parallel_sort.h"


// LDshrinkWriter::LDshrinkWriter(const std::vector<int> &inp,
//                                const LDshrinkWriter::name_vec names,
//                                const bool write_annotations)
//     : snpnames_row(names), snpnames_col(names), p_row(inp.size()), p_col(p_row),
//       write_annotations(write_annotations) {
//   data_store.reserve(((p_row * p_row - p_row) / 2 + p_row) / 5);
//   if (write_annotations) {
//     anno_data_store.reserve(((p_row * p_row - p_row) / 2 + p_row) / 5);
//   }
// }


// LDshrinkWriter::LDshrinkWriter(const std::pair<std::vector<T>,std::vector<T>> &inp,
// 			       const std::pair<LDshrinkWriter::name_vec,LDshrinkWriter::name_vec> names,
// 			       const bool write_annotations):
// 			       snpnames_row(names.first),
// 				 snpnames_col(names.second),
// 				 p_row(inp.first.size()),
// 				 p_col(inp.second.size()),
// 				 write_annotations(write_annotations)
// 			       {
// 				 const size_t tot_p = p_row*p_col;
// 				 data_store.reserve(tot_p/10);
//                                  if (write_annotations) {
//                                    anno_data_store.reserve(tot_p / 10);
//                                  }
// 			       }

// void LDshrinkWriter::write_symm(const std::pair<T,T> idx){
//   data_store.push_back(MatTup<T>(idx));
//  if (write_annotations) {
//   anno_data_store.push_back(MatTup<T>(idx, 0));
//  }
// }
// void LDshrinkWriter::write_nsymm(const std::pair<T,T> idx,const double res,const double annom){
//   data_store.push_back(MatTup<T>(idx,res));
//   if (write_annotations) {
//     anno_data_store.push_back(MatTup<T>(idx, annom));
//   }
// }




// template<typename T, typename F, int N, int Default_a,int Default_b>
// bool operator<(const MatTup<T,F,N,Default_a> &a, const MatTup<T,F,N,Default_b> &b){
//   return(a.m_idx<b.m_idx);
// }


// template<typename T, typename F, int N, int Default_a,int Default_b>
// bool operator==(const MatTup<T,F,N,Default_a> &a, const MatTup<T,F,N,Default_b> &b){
//   return(a.m_idx==b.m_idx);
// }






// void LDshrinkWriter::finalize(){
//   if(write_annotations){
//     const size_t pt=data_store.size();
//     if (pt != anno_data_store.size()) {
//       Rcpp::stop("annotation result size is not equal to LD result size in "
//                  "LDshrinkWriter::finalize!");
//     }
//     tbb::parallel_sort(anno_data_store.begin(),anno_data_store.end());
//     tbb::parallel_sort(data_store.begin(),data_store.end());
//  {
//       for(int i=0; i<pt; i++){
// 	auto ra=data_store[i];
// 	auto rb= anno_data_store[i];
//         if ((ra.row() != rb.row()) || (ra.col() != rb.col())) {
//           Rcpp::Rcerr << "in position: " << i << std::endl;
//           Rcpp::Rcerr << "data_store:" << ra.row() << "," << ra.col() << ","
//                       << ra.value() << std::endl;
//           Rcpp::Rcerr << "anno_data_store:" << rb.row() << "," << rb.col()
//                       << "," << rb.value() << std::endl;
//           Rcpp::stop(
//               "unequal elements in	anno_store and anno_data_store");
//         }
//       }
//     }
//   }
// }

//  Rcpp::DataFrame LDshrinkWriter::toAnnotationDataFrame(const bool use_rownames)const {
//    using namespace Rcpp;
//    const size_t num_trip = data_store.size();
//    const bool stringsAsFactors = false;
//    if(!use_rownames){
//      Vector<INTSXP> rowid(num_trip);
//      Vector<INTSXP> colid(num_trip);
//      NumericVector rvec(num_trip);
//      NumericVector anno(num_trip);
//      for (size_t i = 0; i < num_trip; i++) {
//        auto trip = data_store[i];
//        auto atrip = anno_data_store[i];
//        rowid(i) = trip.row();
//        colid(i) = trip.col();
//        rvec(i) = trip.value();
//        anno(i) = atrip.value();
//      }

//    return (DataFrame::create(_["rowsnp"] = rowid, _["colsnp"] = colid,
// 			     _["r"] = rvec, _["anno"]= anno,
// 			     _["stringsAsFactors"] = stringsAsFactors));
//  }else{
//     Vector<STRSXP> rowid(num_trip);
//     Vector<STRSXP> colid(num_trip);
//     NumericVector rvec(num_trip);
//     NumericVector anno(num_trip);
//     for (size_t i = 0; i < num_trip; i++) {
//       auto trip = data_store[i];
//       auto atrip = anno_data_store[i];
//       rowid(i) = snpnames_row(trip.row());
//       colid(i) = snpnames_col(trip.col());
//       rvec(i) = trip.value();
//       anno(i) = atrip.value();
//     }
//     return (DataFrame::create(_["rowsnp"] = rowid, _["colsnp"] = colid,
//                               _["r"] = rvec, _["anno"] = anno,
//                               _["stringsAsFactors"] = stringsAsFactors));
//    }
//  }

// Rcpp::DataFrame LDshrinkWriter::toDataFrame(const bool use_rownames)const {
//   using namespace Rcpp;

//   if(write_rownnotations){
//     return(toAnnotationDataFrame(use_rownames));
//   }

//   const size_t num_trip = data_store.size();
//   const bool stringsAsFactors = false;

//   if(!use_rownames){
//     Vector<INTSXP> rowid(num_trip);
//     Vector<INTSXP> colid(num_trip);
//     NumericVector rvec(num_trip);
//     for (size_t i = 0; i < num_trip; i++) {
//       auto trip = data_store[i];
//       rowid(i) = trip.row();
//       colid(i) = trip.col();
//       rvec(i) = trip.value();
//     }
//     return (DataFrame::create(_["rowsnp"] = rowid, _["colsnp"] = colid,
// 			      _["r"] = rvec,
// 			      _["stringsAsFactors"] = stringsAsFactors));
//   }else{
//     Vector<STRSXP> rowid(num_trip);
//     Vector<STRSXP> colid(num_trip);
//     NumericVector rvec(num_trip);
//     for (size_t i = 0; i < num_trip; i++) {
//       auto trip = data_store[i];
//       rowid(i) = snpnames_row(trip.row());
//       colid(i) = snpnames_col(trip.col());
//       rvec(i) = trip.value();
//     }
//     return (DataFrame::create(_["rowsnp"] = rowid, _["colsnp"] = colid,
//                               _["r"] = rvec,
//                               _["stringsAsFactors"] = stringsAsFactors));
//   }
// }


// Eigen::SparseMatrix<double> LDshrinkWriter::sparseMatrix() const {
//   static_assert(std::is_integral<T>::value,"Integral index type required for conversion to sparseMatrix");
//   const size_t num_elem = data_store.size();
//   Eigen::SparseMatrix<double> object(p_row, p_col);
//   // Rcpp::Rcerr << "About to create sparse matrix of size: ";
//   // Rcpp::Rcerr << num_elem << std::endl;
//   object.setFromTriplets(data_store.begin(), data_store.end());
//   //    Rcpp::Rcerr << "Sparse Matrix created" << std::endl;
//   return (object);
// }

// SEXP LDshrinkWriter::dsCMatrix() const {
//   auto object = this->sparseMatrix();
//   using namespace Rcpp;
//   const int nnz = object.nonZeros();
//   //  Rcpp::Rcerr << "Number of entries: " << nnz << std::endl;
//   S4 ans("dsCMatrix");
//   ans.slot("Dim") = Dimension(object.rows(), object.cols());
//   ans.slot("i") =
//       IntegerVector(object.innerIndexPtr(), object.innerIndexPtr() + nnz);
//   ans.slot("p") = IntegerVector(
//       object.outerIndexPtr(), object.outerIndexPtr() + object.outerSize() + 1);
//   ans.slot("x") = NumericVector(object.valuePtr(), object.valuePtr() + nnz);
//   ans.slot("Dimnames") = List::create(snpnames_row, snpnames_col);
//   return (ans);
// }

// SEXP LDshrinkWriter::toType(const std::string output_type)const{
//   if(output_type=="dsCMatrix"){
//     return(dsCMatrix());
//   }
//   return(toDataFrame());
// }



// // Rcpp::DataFrame LDshrinkAnnoWriter::toDataFrame()const {
// //   const size_t num_trip = data_store.size();
// //   Rcpp::Vector<cpp2r<T>::data_t> rowid(num_trip);
// //   Rcpp::Vector<cpp2r<T>::data_t> colid(num_trip);
// //   Rcpp::NumericVector rvec(num_trip);
// //   Rcpp::NumericVector avec(num_trip);

//   for (size_t i = 0; i < num_trip; i++) {
//     auto trip = data_store[i];
//     rowid(i) = trip.row();
//     colid(i) = trip.col();
//     rvec(i) = trip.value();
//     avec(i) = trip.anno_value();
//   }
//   using namespace Rcpp;
//   const bool stringsAsFactors = false;
//   return (Rcpp::DataFrame::create(_["rowsnp"] = rowid, _["colsnp"] = colid,
//                                   _["r"] = rvec, _["anno"] = avec,
//                                   _["stringsAsFactors"] = stringsAsFactors));
// }
