#' Calculate LDshrink adjusted LD matrix
#'
#' @param haplo_panel `n` by `p` genotype or haplotype numeric matrix from a contiguous region of the genome ()
#' @param map_data vector of cumulative genetic map values.  Must be 
#'
#' @return LD matrix
#' @export
#'
#' @examples
LDshrink <- function(genotype_panel,map_data,m=85,Ne=11490.672741,cutoff=1e-3,isGeno=NA,cov_2_cor=TRUE,na.rm=TRUE){
  if(is.na(isGeno)){
    isGeno <- max(genotype_panel,na.rm = na.rm)>1
  }
  stopifnot(!is.na(isGeno))
  Genomult <- ifelse(isGeno,0.5,1)
  # haplo_panel <- scale(haplo_panel,center=T,scale=F)
  if(!na.rm){
      return(fastLDshrink(genotype_panel,map_data,m,Ne,cutoff,isGeno,cov_2_cor))
  }else{
    if(cov_2_cor){
      return(stats::cov2cor(shrinkCov(stats::cov(genotype_panel,use=ifelse(na.rm,"complete.obs","all.obs"))*Genomult,map_data,m,Ne,cutoff)))
    }else{
      return(shrinkCov(stats::cov(genotype_panel,use=ifelse(na.rm,"complete.obs","all.obs"))*Genomult,map_data,m,Ne,cutoff))
    }
  }
}
