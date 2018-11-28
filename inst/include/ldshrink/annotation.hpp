#pragma once

#include <RcppEigen.h>


class ConstantGeneticMap{
  const bool strict;
  using	idmap=std::map<int,double>;
  const idmap genmap;
public:
  ConstantGeneticMap(const Rcpp::IntegerVector &pos_vec,
                     const Rcpp::NumericVector &map_vec,
		     const bool strict_);
  double interpolate(const int pos)const;
private:
  double linear_interp(const idmap::const_iterator prev_p,const idmap::const_iterator next_p, const int pos) const;
  double interpolate_prev(const int pos) const;
  double interpolate_post(const int pos) const;

};



// Rcpp::XPtr<ConstantGeneticMap>
// new_constantmap(const Rcpp::NumericVector map_data,
//                 const Rcpp::IntegerVector map_pos);


class DistAnnoVec{
  const std::vector<double> annot;
public:
  DistAnnoVec(const Rcpp::NumericVector genmap) : annot(Rcpp::as<std::vector<double> >(genmap)) {}
  double get(const std::array<int,2> idx) const {
    return (std::fabs(annot[idx[0]] - annot[idx[1]]));
  }
};

class DistAnno {
  using T = int;
  using snpmap = std::unordered_map<T, double>;
  const snpmap map_data;

public:
  DistAnno(const Rcpp::NumericVector genmap)
    : map_data([&]() {
		 Rcpp::RObject data_attr = genmap.attr("index");
		 std::vector<T> ret;
		 const size_t p = genmap.size();
		 snpmap retmap;
		 retmap.reserve(p);
		 if (data_attr.isNULL()) {
		   data_attr = genmap.attr("names");
		   if (data_attr.isNULL()) {
		     Rcpp::stop("annotation must be a named (character) vector, or "
				"names must be provided (there must be an "
				"attributes 'names' or 'index'");
		   }
		 }
		 ret = Rcpp::as<std::vector<T>>(data_attr);
		 double omap = genmap[0];
		 for (size_t i = 0; i < p; i++) {
		   double tmap = genmap(i);
		   T idxn = ret[i];
		   if (i > 0) {
		     if (tmap <= omap) {
		       Rcpp::Rcerr << "In position :" << i + 1
				   << " of annotation:" << std::endl;
		       Rcpp::stop("annotation must be (strictly) sorted");
		     }
		   }
		   retmap.insert({idxn, tmap});
		 }
		 return (std::move(retmap));
	       }()) {}
  double get(const std::pair<T, T> idx) const {
    return (std::fabs(map_data.at(idx.first) - map_data.at(idx.second)));
  }
};







// template<typename T>
// class ConstDistAnno{
//   using	snpmap=std::unordered_map<T,double>;
//   const Rcpp::XPtr<ConstantGeneticMap> map_ctor;
//   const snpmap map_data;

// public:
//   DistAnno(const Rcpp::IntegerVector pos,
// 	   Rcpp::XPtr<ConstantGeneticMap> ctor):map_ctor(ctor),
// 						map_data([&](){
// 							   Rcpp::RObject data_attr=genmap.attr("index");
// 							   std::vector<T> ret;
// 							   const size_t p=genmap.size();
// 							   snpmap retmap;
// 							   retmap.reserve(p);
// 							   if(data_attr.isNULL()){
// 							     data_attr=pos.attr("names");
// 							     if(data_attr.isNULL()){
// 							       Rcpp::stop("annotation must be a named (character) vector, or names must be provided (there must be an attributes 'names' or 'index'");
// 							     }
// 							   }
// 							   ret =Rcpp::as<std::vector<T> >(data_attr);
// 							   double omap= genmap[0];
// 							   for(size_t i=0; i<p; i++){
// 							     double tmap=ctor.interpolate(pos[i]);
// 							     T idxn=ret[i];
// 							     if(i>0){
// 							       if(tmap<=omap){
// 								 Rcpp::Rcerr<<"In position :"<< i+1 <<" of annotation:"<<std::endl;
// 								 Rcpp::stop("annotation must be (strictly) sorted");
// 							       }
// 							     }
// 							     retmap.insert({idxn,tmap});
// 							   }
// 							   return(std::move(retmap));
// 							 }()){
//   }
//   double get(const std::pair<T,T> idx)const{
//     return(std::fabs(map_data.at(idx.first)-map_data.at(idx.second)));
//   }
// };
