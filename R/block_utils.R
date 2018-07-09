#' Assign SNPs to LD blocks
#'
#' @param snp_df dataframe of snp coordinates, must contain columns named `chr` and `pos`.
#' @param break_df dataframe of LD blocks, if `NULL`, use precomputed (EUR) ldshrink LD blocks, must contain columns named `chr` `start`, `stop`, and `region_id`
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#'
#' @return modified `snp_df` dataframe with additional column `region_id` mapping snp to LD block
#' @export
#'
#' @examples
assign_snp_block <- function(snp_df,break_df=NULL,assign_all=T){
  stopifnot(dplyr::group_by(snp_df,chr) %>% dplyr::summarise(all_sort=all(!is.unsorted(pos,strictly = T))) %>% dplyr::summarise(as=all(all_sort)) %>% dplyr::pull(as),
            min(snp_df$pos)>0,is.integer(snp_df$chr),is.integer(snp_df$pos))
  if(is.null(break_df)){
      liftover_allf <- system.file("fourier_ls-all.bed.gz",package="LDshrink")
      break_df <- readr::read_delim(liftover_allf,delim="\t",trim_ws = T,col_names=c("chr","start","stop"))
      break_df <- dplyr::mutate(break_df,chr=as.integer(gsub("chr","",chr))) %>% dplyr::mutate(region_id=1:n())
  }else{
   stopifnot(!is.null(break_df$chr),
             !is.null(break_df$start),
             !is.null(break_df$stop),
             !is.null(break_df$region_id),
             is.integer(break_df$chr),
             is.integer(break_df$start),
             is.integer(break_df$stop),
     dplyr::group_by(break_df,chr) %>% dplyr::summarise(all_sort=all((!is.unsorted(start,strictly = T)) & (!is.unsorted(stop,strictly = T)))) %>% dplyr::summarise(as=all(all_sort)) %>% dplyr::pull(as),
             all(break_df$start<break_df$stop)
             )
  }
  
  return(dplyr::mutate(snp_df,region_id=set_ld_region(break_df,snp_df,assign_all = assign_all)))
}






#' Assign or Interpolate Genetic Map Values
#'
#' @param snp_df 
#' @param map_df dataframe of LD blocks, if `NULL`, use precomputed (EUR) ldshrink LD blocks
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#'
#' @return
#' @export
#'
#' @examples
assign_map <- function(snp_df,map_df){
  u_chr <- dplyr::distinct(snp_df,chr)
  snp_dfl <- split(snp_df,snp_df$chr)
  map_dfl <- dplyr::semi_join(map_df,u_chr,by="chr") %>% split(.$chr)
  stopifnot(all(names(map_dfl)==names(snp_dfl)))
  retdf <- purrr::map2_df(map_dfl,snp_dfl,~dplyr::mutate(.y,map=interpolate_map(.x$map,.x$pos,.y$pos)))
  return(retdf)
}


#' Download 1000 genomes OMNI genetic map data
#'
#' @param pop population, see `Details`
#' @param destination local pathname for destination
#'
#' @details `pop` can be one of `"ACB" "ASW" "CDX" "CEU" "CHB" "CHS" "CLM" "FIN" "GBR" "GIH" "IBS" "JPT" "KHV" "LWK" "MKK" "MXL" "PEL" "PUR" "REA" "TSI" "YRI"`
#' @return pathname of downloaded file
#' @export
#'
#' @examples 
download_omni_map <- function(pop="CEU", destination_dir =tempdir()){
  destination <- tempfile()
  base_url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/"
  download.file(paste0(base_url,pop,"_omni_recombination_20130507.tar"),destfile = destination)
  untar(destination,exdir=destination_dir)
  result_files <- dir(file.path(destination_dir,pop),full.names = T)
  return()
}





LDshrink_evd <- function(panel,map=NULL,m=85,
                            Ne=11490.672741,
                            cutoff=1e-3,
                         useLDshrink=T,na.rm=F){
#    stopifnot(ncol(panel)==len
    isGeno <- max(panel,na.rm = na.rm)>1
    if(useLDshrink){
        S <- LDshrink(genotype_panel = panel,
                      map_data = map,
                      m = m,
                      Ne = Ne,
                      cutoff = cutoff,
                      cov_2_cor = T,
                      na.rm = na.rm)
    }else{
      S <- stats::cor(panel,use = "complete.obs")
    }
    L2 <- colSums(S^2)-1
    evdR <- eigen(S)
    return(list(R=S,L2=L2,D=evdR$values,Q=evdR$vectors))
}
LDshrink_df <- function(panel, map,
                        snp_id,
                        m=85,
                        Ne=11490.672741,
                        cutoff=1e-3,
                        r2_cutoff=0.01,
                        useLDshrink=T,
                        na.rm=F, progress=FALSE){

  if(useLDshrink){
      return(ld2df_p(scaled_data=scale(panel,center=T,scale=F),
                     mapd=map,
                     m = m,
                     Ne = Ne,
                     cutoff = cutoff,
                     rsid = snp_id,
                     r2cutoff = r2_cutoff,progress=progress))

  }else{
    return(ld2df(ldmat=stats::cor(panel,use = "complete.obs"),
          rsid=snp_id,
          r2cutoff=r2_cutoff))
  }
}




## flip_allele_exp <- function(allele_a,allele_b){
##   utf_i <- Vectorize(utf8ToInt)
##   data_df <- tibble::data_frame(allele_a=allele_a,allele_b=allele_b)
##   gwas_snp_df <- data_df %>% 
##     tidyr::separate(allele_a,c("ref_a","alt_a")) %>% 
##     tidyr::separate(allele_b,c("ref_b","alt_b")) %>%
##     dplyr::mutate(
##       ref_a=utf_i(tolower(ref_a)),
##       alt_a=utf_i(tolower(alt_a)),
##       ref_b=utf_i(tolower(ref_b)),
##       alt_b=utf_i(tolower(alt_b))
##     )
 
##   return(flip_allele(gwas_snp_df$ref_a,
##                      gwas_snp_df$alt_a,
##                      gwas_snp_df$ref_b,
##                      gwas_snp_df$alt_b))
## }


