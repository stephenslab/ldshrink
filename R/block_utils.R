#' Assign SNPs to LD blocks
#'
#' @param snp_df 
#' @param break_df dataframe of LD blocks, if `NULL`, use precomputed (EUR) ldshrink LD blocks
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#'
#' @return
#' @export
#'
#' @examples
assign_snp_block <- function(snp_df,break_df=NULL,assign_all=T){
  if(is.null(break_df)){
    data("break_df")
  }
  return(dplyr::mutate(snp_df,region_id=set_ld_region(break_df,snp_df,assign_all = assign_all)))
}


assign_map <- function(snp_df,map_df){
  u_chr <- dplyr::distinct(snp_df,chr)
  snp_dfl <- split(snp_df,snp_df$chr)
  map_dfl <- dplyr::semi_join(map_df,u_chr,by="chr") %>% split(.$chr)
  stopifnot(all(names(map_dfl)==names(snp_dfl)))
  retdf <- purrr::map2_df(map_dfl,snp_dfl,~dplyr::mutate(.y,map=interpolate_map(.x$map,.x$pos,.y$pos)))
  return(retdf)
}



LDshrink_evd <- function(panel,map=NULL,m=85,
                            Ne=11490.672741,
                            cutoff=1e-3,
                         useLDshrink=T,na.rm=F){
#    stopifnot(ncol(panel)==len
    isGeno <- max(panel,na.rm = na.rm)>1
    if(useLDshrink){
        S <- LDshrink(haplo_panel = panel,
                      map_data = map,
                      m = m,
                      Ne = Ne,
                      cutoff = cutoff,
                      cov_2_cor = T,
                      na.rm = na.rm)
    }else{
      S <- cor(panel,use = "complete.obs")
    }
    L2 <- colSums(S^2)-1
    evdR <- eigen(S)
    return(list(R=S,L2=L2,D=evdR$values,Q=evdR$vectors))
}
    






