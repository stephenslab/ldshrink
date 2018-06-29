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
#include <testthat.h>
#include "LDshrink.h"


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Shrinkage functor") {
  const	int p=3;
  const double m=85;
  const double Ne=11490.672741;
  const double cutoff=0.001;

  Eigen::VectorXd map_data(p);

  map_data<<0.317759,0.394240,1.067599;
  Eigen::Map<const Eigen::VectorXd> map_map(map_data.data(),p);

  Eigen::MatrixXd true_map_matrix(p,p);
  Eigen::MatrixXd mapfmat = makeMapDiff(map_map,m,Ne,cutoff);
  for(int i=0;i<p;i++){
    for(int j=0; j<p;j++){
      const double map_dist=std::fabs(map_data(j)-map_data(i));
      double rho=4*Ne*map_dist/100;
      rho=-rho/(2*m);
      rho=std::exp(rho);
      true_map_matrix(i,j)=rho<cutoff ? 0 : rho;
      expect_true(mapfmat(i,j)==true_map_matrix(i,j));

    }


  }
  expect_true(mapfmat.isApprox(true_map_matrix));
}


//   // The format for specifying tests is similar to that of
//   // testthat's R functions. Use 'test_that()' to define a
//   // unit test, and use 'expect_true()' and 'expect_false()'
//   // to test the desired conditions.
//   // Rcpp::IntegerVector gwas_ref= {0,0,1,1,2,2};
//   // Rcpp::IntegerVector gwas_alt= {1,1,2,2,1,1};

//   // Rcpp::IntegerVector ld_ref= {0,0,1,1,2,2};
//   // Rcpp::IntegerVector ld_alt= {1,1,2,2,1,1};




//   // test_that("flip_allele doesn't flip any when they're the same") {
//   //   expect_true(flip_allele(gwas_ref,gwas_alt,ld_ref,ld_alt)==Rcpp::LogicalVector::create(true,true,true

//   // }

// }
