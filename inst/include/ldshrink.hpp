#pragma once
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
//#include "xtensor/xmath.hpp"          // xtensor import for the C++ universal functions
//#include "xtensor-r/rarray.hpp"       // R bindings

//[[Rcpp]]
#include<RcppCommon.h>

namespace Rcpp{
  namespace traits{
    //template<> std::vector<std::string_view> as(SEXP);
    SEXP wrap(const std::vector<std::string_view> & obj);
    template <> class Exporter< std::vector<std::string_view> >;
    template <> class Exporter<std::string_view>;
  }
}




#include <RcppEigen.h>


#include "ldshrink/LD.hpp"
#include "ldshrink/misc.hpp"
#include "ldshrink/ld_writer.hpp"
#include "ldshrink/genotype.hpp"
#include "ldshrink/annotation.hpp"
#include "ldshrink/cov.hpp"
#include "ldshrink/parameter_list.hpp"
