#include <RcppEigen.h>
#include<algorithm>
#include <functional>
#include <tuple>
#include "ldshrink.hpp"



//[[Rcpp::export]]
SEXP round_trip_skyline_t(Rcpp::NumericMatrix r_mat,const std::string to_t){

  std::array<Rcpp::NumericMatrix, 1> r_array{r_mat};

  //  Rcpp::Rcout<<"r_mat:"<<r_mat<<std::endl;
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  Skyline_data_store<1> test_store(r_array,true);
  test_store.finalize();
  return (test_store.toType(to_t, rownames(r_mat), colnames(r_mat)));
}
