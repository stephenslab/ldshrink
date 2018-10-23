#'




#' Calculate ldshrink adjusted LD matrix
#'
#' @param genotype_panel `n` by `p` genotype or haplotype numeric matrix from a contiguous region of the genome ()
#' @param map_data vector of cumulative genetic map values.  Must be sorted
#'
#' @return LD matrix
#' @export
#'
ldshrink <- function(genotype_panel, map_data, m=85, Ne=11490.672741, cutoff=1e-3, isGeno=NA, cov_2_cor=TRUE, na.rm=TRUE, process_map=TRUE){
  if(is.na(isGeno)){
    isGeno <- max(genotype_panel, na.rm = na.rm)>1
  }
  stopifnot(!is.na(isGeno))
  Genomult <- ifelse(isGeno, 0.5, 1)
  # haplo_panel <- scale(haplo_panel, center=T, scale=F)
  if(!na.rm){
      return(fastldshrink(genotype_panel, map_data, m, Ne, cutoff, isGeno, cov_2_cor))
  }else{
    if(cov_2_cor){
      return(stats::cov2cor(shrinkCov(stats::cov(genotype_panel, use=ifelse(na.rm, "complete.obs", "all.obs"))*Genomult, map_data, m, Ne, cutoff)))
    }else{
      return(shrinkCov(stats::cov(genotype_panel, use=ifelse(na.rm, "complete.obs", "all.obs"))*Genomult, map_data, m, Ne, cutoff))
    }
  }
}



ldshrink_evd <- function(panel, map=NULL, m=85,
                            Ne=11490.672741, 
                            cutoff=1e-3, 
                         useldshrink=T, na.rm=F){
#    stopifnot(ncol(panel)==len
    isGeno <- max(panel,na.rm = na.rm)>1
    if(useldshrink){
        S <- ldshrink(genotype_panel = panel,
                      map_data = map,
                      m = m,
                      Ne = Ne,
                      cutoff = cutoff,
                      cov_2_cor = T,
                      na.rm = na.rm)
    }else{
      S <- stats::cor(panel, use = "complete.obs")
    }
    N <- nrow(panel)
    L2 <- colSums(S^2)-1
    L2 <- L2-(1-L2)/(N-2)
    evdR <- eigen(S)
    return(list(R=S, L2=L2, D=evdR$values, Q=evdR$vectors))
}
