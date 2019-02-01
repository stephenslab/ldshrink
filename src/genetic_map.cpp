#include <map>
#include <ldshrink.hpp>
#include <Rcpp.h>
#include <progress.hpp>
#include <ldshrink/annotation.hpp>

double ConstantGeneticMap::interpolate_post(const int pos)const {
  const std::string error_mess="position"+std::to_string(pos)+" is after final position";
  if (strict) {
    Rcpp::stop(error_mess);
  }
  //  Rcpp::Rcerr << error_mess << std::endl;
  const auto it_end = std::prev(genmap.end());
  const auto it_prev = std::prev(it_end);
  return (linear_interp(it_prev, it_end, pos));
}

double ConstantGeneticMap::interpolate_prev(const int pos)const{
  const std::string error_mess="position"+std::to_string(pos)+" is before first position";
  if (strict) {
    Rcpp::stop(error_mess);
  }
  //  Rcpp::Rcerr << error_mess << std::endl;
  const auto gb=genmap.begin();
  if (gb->first == pos) {
    return (gb->second);
  }
  const auto gbp = std::next(gb);
  return (linear_interp(gb, gbp, pos));
}

double ConstantGeneticMap::linear_interp(const idmap::const_iterator prev_it,
                                          const idmap::const_iterator next_it,
                                          const int pos) const {

  const int prev_mappos = prev_it->first;
  const int mappos = next_it->first;
  const double prev_map = prev_it->second;
  const double cur_map = next_it->second;
  double frac = static_cast<double>(pos - prev_mappos) /
    static_cast<double>(mappos - prev_mappos);
  if(strict){
    if (frac < 0 || frac > 1) {
      Rcpp::Rcerr << "in position: " << pos << std::endl;
      Rcpp::Rcerr << prev_it->first << "," << prev_it->second << std::endl;
      Rcpp::Rcerr << next_it->first << "," << next_it->second << std::endl;
      Rcpp::Rcerr << "frac:" << frac << std::endl;
      Rcpp::stop("frac must be positive (and less than 1)");
    }
  }
  return (prev_map + frac * (cur_map - prev_map));
}

double ConstantGeneticMap::interpolate(const int pos)const {
  //Return first value that is greater than or equal to pos
  auto n_it = genmap.lower_bound(pos);
  if(n_it == genmap.end()){
    return (interpolate_post(pos));
  }
  if(n_it == genmap.begin()){

    return (interpolate_prev(pos));
  }
  if (n_it->first == pos) {
    return (n_it->second);
  }
  const auto prev_it = std::prev(n_it);
  return (linear_interp(prev_it, n_it, pos));
}

ConstantGeneticMap::ConstantGeneticMap(const Rcpp::IntegerVector &pos_vec,
                                       const Rcpp::NumericVector &map_vec,
                                       const bool strict_)
  :       strict(strict_),
	  genmap([](const Rcpp::IntegerVector &pos_vec,
		    const Rcpp::NumericVector &map_vec,
		    const bool strict) {
		   auto pos_begin = pos_vec.begin();
		   auto pos_end = pos_vec.end();
		   auto map_begin = map_vec.begin();
		   auto map_end = map_vec.end();
		   const size_t p = pos_end - pos_begin;
		   if (p != (map_end - map_begin)) {
		     Rcpp::stop("position and genetic map must be the same size");
		   }
		   if (p < 2) {
		     Rcpp::stop("There must be at least 2 elements in the genetic map in "
				"order to interpolate!");
		   }
		   std::map<int, double> retmap;
		   auto rb = retmap.begin();
		   auto pb = pos_begin;
		   auto mb = map_begin;
		   double orb=*mb;
		   rb = retmap.insert(rb, std::make_pair(*(pb), *(mb)));
		   for (size_t i = 1; i < p; i++) {
		     pb++;
		     mb++;
		     if (orb >= *mb) {
		       Rcpp::Rcerr << "Genetic map is not strictly sorted at position: "
				   << i << " (" << orb << ">=" << *mb << ")" << std::endl;
		       if (strict || orb > *mb) {
			 Rcpp::stop("Genetic map must be strictly sorted ");
		       }
		     } else {
		       rb = retmap.insert(rb, std::make_pair(*(pb), *(mb)));
		     }
		     orb = *mb;
		   }
		   return (std::move(retmap));
		 }(pos_vec, map_vec, strict_)) {}


//' Linear interpolation of genetic map values
// '@param map  is a length `p` vector of cumulative genetic map values. `map` must be _strictly_ _sorted_
// '@param map_pos  is a length `p` vector of genome coordinates corresponding to the reference genetic map. `map_pos` must be _strictly_ _sorted_
// '@param target_pos is a vector of coordinates to interpolate
// '@param strict a boolean indicating whether to strictly interpolate
// '@param progress a boolean indicating whether to indicate progress with a progress bar
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector interpolate_genetic_map(const Rcpp::NumericVector &map,
					    const Rcpp::IntegerVector map_pos,
					    const Rcpp::IntegerVector target_pos,
					    const bool strict = true,
					    const bool progress = false) {

  const size_t p = target_pos.size();
  ConstantGeneticMap ref_map(map_pos, map, strict);
  Rcpp::NumericVector ret(p);
  Progress prog_bar(p, progress);
  std::transform(target_pos.begin(), target_pos.end(), ret.begin(),
                 [&](const int i) {
		   prog_bar.increment();
		   return (ref_map.interpolate(i)); });
  return (ret);
}

