/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.

#include "ldshrink.hpp"

#include <Rcpp.h>
#include <testthat.h>
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.
Rcpp::NumericMatrix rand_mat(const size_t a, const size_t b,const bool symmetric=false) {
  using namespace Rcpp;
  RNGScope rngScope;

  NumericMatrix retmat(a,b);
  if (symmetric) {
    for(int i=0; i<a; i++){
      for (int j = i; j < b; j++) {
        retmat(j, i) = R::runif(0, 1);
        retmat(i, j) = retmat(j, i);
      }
    }
  }else{
    for(int i=0; i<a; i++){
      for (int j = 0; j < b; j++) {
        retmat(j, i) = R::runif(0, 1);
      }
    }
  }

  return (retmat);
}

Rcpp::StringVector gen_name(const size_t n){
  Rcpp::StringVector row_id(n);
  for (int i = 0; i < n; i++) {
    row_id(i) = Rcpp::wrap(std::to_string(i));
  }
  return (row_id);
}

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Skyline roundtrip") {

  const size_t p_a=10;
  const size_t p_b=10;


  auto r_mat =rand_mat(p_a,p_b,true);
  auto r_id= gen_name(p_a);
  rownames(r_mat)=r_id;
  colnames(r_mat)=r_id;

  std::array<Rcpp::NumericMatrix, 1> r_array { r_mat };

  //  Rcpp::Rcout<<"r_mat:"<<r_mat<<std::endl;
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  Skyline_data_store<1> test_store(r_array);
  test_that("roundtrip from symmetric matrix") {
    Rcpp::NumericMatrix ret = test_store.toMatrix(r_id,r_id);
    //Rcpp::Rcout<<"ret:"<<ret<<std::endl;
    for(int i=0; i<p_a; i++){
      for (int j = 0; j < p_b; j++) {
        expect_true(ret(i, j) == r_mat(i, j));
      }
    }
  }
}
