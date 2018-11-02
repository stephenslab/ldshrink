#include <map>
#include "ldshrink.hpp"
#include <Rcpp.h>

// [[Rcpp::depends(BH)]]


double ConstantGeneticMap::interpolate(const int pos)const {
  auto n_it = genmap.lower_bound(pos);
  if (n_it == genmap.end() || n_it==genmap.begin()) {
    if(n_it->first==pos){
      return(n_it->second);
    }
    const std::string error_mess="position"+std::to_string(pos)+" is before (or at) first position";
    if (strict) {
      Rcpp::stop(error_mess);
    }
    Rcpp::Rcerr << error_mess << std::endl;
    return (0);
  }
  if(n_it->first==pos){
    return(n_it->second);
  }


  const auto prev_it =	std::prev(n_it);
  if (n_it == genmap.end()) {
    const std::string error_mess =
      "position" + std::to_string(pos) + " is after final position";
    Rcpp::stop(error_mess);
  }
  const int prev_mappos = prev_it->first;
  const int mappos = n_it->first;
  const double prev_map = prev_it->second;
  const double cur_map = n_it->second;


  double frac = static_cast<double>(pos - prev_mappos) /
    static_cast<double>(mappos - prev_mappos);

  if(frac<0 || frac>1){
    Rcpp::Rcerr<<"in position: "<<pos<<std::endl;
    Rcpp::Rcerr<<prev_it->first<<","<<prev_it->second<<std::endl;
    Rcpp::Rcerr<<n_it->first<<","<<n_it->second<<std::endl;
    Rcpp::Rcerr<<"frac:"<<frac<<std::endl;
    Rcpp::stop("frac must be positive (and less than 1)");
  }

  return (prev_map + frac * (cur_map - prev_map));
}

ConstantGeneticMap::ConstantGeneticMap(const Rcpp::IntegerVector &pos_vec,
                                       const Rcpp::NumericVector &map_vec,
                                       const bool strict_)
    : genmap([](const Rcpp::IntegerVector &pos_vec,
                const Rcpp::NumericVector &map_vec) {
        auto pos_begin = pos_vec.begin();
        auto pos_end = pos_vec.end();
        auto map_begin = map_vec.begin();
        auto map_end = map_vec.end();
        const size_t p = pos_end - pos_begin;
        if (p != (map_end - map_begin)) {
          Rcpp::stop("position and genetic map must be the same size");
        }
        std::map<int, double> retmap;
        auto rb = retmap.begin();
        auto pb = pos_begin;
        auto mb = map_begin;
        for (size_t i = 0; i < p; i++) {
          rb = retmap.insert(rb, std::make_pair(*(pb++), *(mb++)));
        }
        return (std::move(retmap));
      }(pos_vec, map_vec)),
      strict(strict_) {}




//[[Rcpp::export]]
Rcpp::NumericVector interpolate_map(const Rcpp::NumericVector &map,const Rcpp::IntegerVector map_pos,const Rcpp::IntegerVector target_pos,const bool progress=false){

    const size_t map_snps = map.size();
    if (map_pos.size() != map_snps) {
      Rcpp::stop("map and map_pos must be the same size!");
    }
    const size_t t_snps = target_pos.size();
    Rcpp::NumericVector retvec(t_snps);
    if (!std::is_sorted(map.begin(), map.end())) {
      Rcpp::stop("Genetic map must be sorted!");
    }
    if (!std::is_sorted(map_pos.begin(), map_pos.end(),
                        [](const int &a, const int &b) { return (a < b); })) {
      Rcpp::stop("Reference physical map must be sorted!");
    }
    if (!std::is_sorted(target_pos.begin(), target_pos.end(),
                        [](const int &a, const int &b) { return (a < b); })) {
      Rcpp::stop("Reference physical map must be sorted!");
    }

    Progress prog_bar(t_snps, progress);
    size_t idx2 = 0;
    size_t idx1 = 0;
    while (idx1 < t_snps) {
      const int pos = target_pos[idx1];
      const int mappos = map_pos[idx2];
      if (pos == mappos) {
        retvec[idx1] = map[idx2];
        idx1++;
        prog_bar.increment();
      } else {
        if (pos < mappos) {
          if (idx2 == 0) {
            retvec[idx1] = map[idx2];
            idx1++;
          } else {
            double prev_map = map[idx2 - 1];
            double prev_mappos = static_cast<double>(map_pos[idx2 - 1]);
            double frac = (static_cast<double>(pos) - prev_mappos) /
                          (static_cast<double>(mappos) - prev_mappos);
            retvec[idx1] = prev_map + frac * (map[idx2] - prev_map);
            idx1++;
            prog_bar.increment();
          }
        } else {
          if (pos > mappos) {
            if (idx2 == (map_snps - 1)) {
              retvec[idx1] = map[idx2];
              idx1++;
              prog_bar.increment();
            } else {
              idx2++;
            }
          } else {
            Rcpp::stop("Something has gone wrong");
          }
        }
      }
    }
  return(retvec);
  }
