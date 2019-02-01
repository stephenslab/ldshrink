#pragma once
#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>

template <typename T>
struct cpp2r{
 static const SEXPTYPE data_t = NILSXP;
};

template <>
struct cpp2r<bool>{
 static const SEXPTYPE data_t = LGLSXP;
};

template<>
struct cpp2r<int>{
   static const SEXPTYPE data_t =INTSXP;
};

template<>
struct cpp2r<double>{
   static const SEXPTYPE data_t =REALSXP;
};

template<>
struct cpp2r<std::string>{
   static const SEXPTYPE data_t =STRSXP;
};




template<int RTYPE> struct r2cpp_t{
  typedef std::false_type type;
};
template<> struct r2cpp_t<INTSXP>{
  typedef int type;
};
template<> struct r2cpp_t<REALSXP>{
  typedef double type;
};
template<> struct r2cpp_t<LGLSXP>{
  typedef bool type;
};
template<> struct r2cpp_t<CHARSXP>{
  typedef std::string type;
};
#if __cplusplus > 201402L
template<> struct r2cpp_t<STRSXP>{
  typedef std::string_view type;
};
#else
template<> struct r2cpp_t<STRSXP>{
  typedef std::string type;
};
#endif

template<typename T>
inline T value_or(Rcpp::List dat,const std::string name,const T alt){
  auto myn = Rcpp::as<Rcpp::StringVector>(dat.names());
  auto fn = std::find(myn.begin(),myn.end(),name);
  if (fn == myn.end()) {
    return (alt);
  }
  return (Rcpp::as<T>(dat[name]));
}
