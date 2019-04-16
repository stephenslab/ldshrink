#include "ldshrink.hpp"
#include <algorithm>
#include <functional>

#include <progress.hpp>
#include <progress_bar.hpp>
#include <tuple>
//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]
struct SNPpos{
  std::pair<int,int> c_p;
  SNPpos(const int chrom,const int pos):c_p(std::make_pair(chrom,pos)){};
};

struct GRange{
  std::pair<SNPpos,SNPpos> range;
  GRange(const int chrom,const int start,const int end):range(std::make_pair(SNPpos(chrom,start),SNPpos(chrom,end))){
    if(end<start){
      Rcpp::Rcerr<<"For range: "<<chrom<<":"<<start<<"-"<<end<<std::endl;
      Rcpp::stop("Invalid range! end<start");
    }};
};
bool operator<(const SNPpos p, const GRange g){
  return(p.c_p<g.range.first.c_p);
}
bool operator<=(const SNPpos p, const GRange g){
  return(p.c_p<=g.range.first.c_p);
}
bool operator>(const SNPpos p, const GRange g){
  return(p.c_p>g.range.second.c_p);
}
bool operator>=(const SNPpos p, const GRange g){
  return(p.c_p>=g.range.second.c_p);
}
bool operator==(const SNPpos p, const GRange g){
  return((p.c_p>=g.range.first.c_p) && (p.c_p<=g.range.second.c_p));
}




//[[Rcpp::export]]
bool sorted_snp_df(const Rcpp::DataFrame &snp_info){
  const Rcpp::IntegerVector chr= snp_info["chr"];
  const Rcpp::IntegerVector pos= snp_info["pos"];

  const size_t p=chr.size();
    //std::vector<std::pair<int,int>> r(p);
  std::pair<int,int> old_p{chr[0],pos[0]};
  std::pair<int,int> new_p;
  for(int i=1; i<p;i++){
    new_p={chr[i],pos[i]};
    if(new_p<old_p){
      return(false);
    }
    old_p=new_p;
  }
  return(true);
}



//[[Rcpp::export]]
Rcpp::IntegerVector set_ld_region(const Rcpp::IntegerVector ld_chr,
                                  const Rcpp::IntegerVector ld_start,
                                  const Rcpp::IntegerVector ld_stop,
                                  const Rcpp::IntegerVector ld_region_id,
                                  const Rcpp::IntegerVector chr,
                                  const Rcpp::IntegerVector pos,
                                  const bool assign_all=true){

  const size_t ld_size=ld_chr.size();

  if(!std::is_sorted(ld_chr.begin(),ld_chr.end())){
    Rcpp::stop("break regions must be sorted by chromosome!");
  }

  const size_t snp_size = chr.size();
  Rcpp::IntegerVector ret_region(snp_size);
  if(!std::is_sorted(chr.begin(),chr.end())){
    Rcpp::stop("SNPs must be sorted by chromosome!");
  }
  size_t idx1=0;
  size_t idx2=0;
  GRange gr(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
  while(idx1<snp_size){
    
    //    auto range_start= std::make_tuple(ld_chr[idx2],ld_start[idx2]);
    //    auto range_stop= std::make_tuple(ld_chr[idx2],ld_stop[idx2]);
    const SNPpos s_pos(chr[idx1],pos[idx1]);
    if(s_pos==gr){
      //point in interval
      ret_region[idx1]=ld_region_id[idx2];
      idx1++;
    }else{
      if(s_pos<gr){
        // point before interval
        if(assign_all){
          Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
          Rcpp::Rcerr<<"Can't map to "<<ld_chr[idx2]<<":"<<ld_start[idx2]<<"-"<<ld_stop[idx2]<<" (ahead)"<<std::endl;
          if(idx2>0){
            Rcpp::Rcerr<<"Can't map to "<<ld_chr[idx2-1]<<":"<<ld_start[idx2-1]<<"-"<<ld_stop[idx2-1]<<" (behind)"<<std::endl;
          }
          Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
        }else{
          ret_region[idx1]=NA_INTEGER;
        }
        idx1++;
      }
      if(s_pos>gr){
        //point after interval
        if((idx2+1)<ld_size){
          idx2++;
          gr=GRange(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
        }else{
          if(assign_all){
            Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
            Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
          }else{
            ret_region[idx1]=NA_INTEGER;
          }
          idx1++;
        }
      }
    }
  }
  return(ret_region);
}
