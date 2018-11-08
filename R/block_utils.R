#' Assign SNPs to LD blocks
#'
#' @param snp_df dataframe of snp coordinates, must contain columns named `chr` and `pos`.
#' @param break_df dataframe of LD blocks, if `NULL`, use precomputed (EUR) ldshrink LD blocks, must contain columns named `chr` `start`, `stop`, and `region_id`
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#'
#' @return modified `snp_df` dataframe with additional column `region_id` mapping snp to LD block
#' @export
assign_snp_block <- function(snp_df, break_df=NULL, assign_all=T){
  stopifnot(dplyr::group_by(snp_df, chr) %>% dplyr::summarise(all_sort=all(!is.unsorted(pos, strictly = T))) %>% dplyr::summarise(as=all(all_sort)) %>% dplyr::pull(as), 
            min(snp_df$pos)>0, is.integer(snp_df$chr), is.integer(snp_df$pos))
  return(dplyr::mutate(snp_df, region_id=assign_region(break_df, snp_df, assign_all = assign_all)))
}


assign_region <- function(break_coord, snp_coord, assign_all){
  UseMethod("assign_region")
}

assign_region.data.frame <- function(break_coord, snp_coord, assign_all){
  stopifnot(!is.null(break_df$chr), 
            !is.null(break_df$start), 
            !is.null(break_df$stop),
            !is.null(break_df$region_id),
            is.integer(break_df$chr),
            is.integer(break_df$start),
            is.integer(break_df$stop),
            dplyr::group_by(break_df, chr) %>%
              dplyr::summarise(all_sort=all((!is.unsorted(start, strictly = T)) & (!is.unsorted(stop, strictly = T)))) %>%
              dplyr::summarise(as=all(all_sort)) %>% dplyr::pull(as),
            all(break_df$start<break_df$stop)
  )
  return(set_ld_region(break_coord, snp_coord, assign_all=assign_all))
}



assign_region.function <- function(break_coord, snp_coord, assign_all){
  break_df <- break_coord()
  stopifnot(!is.null(break_df$chr),
            !is.null(break_df$start),
            !is.null(break_df$stop),
            !is.null(break_df$region_id),
            is.integer(break_df$chr),
            is.integer(break_df$start),
            is.integer(break_df$stop),
            dplyr::group_by(break_df, chr) %>% 
              dplyr::summarise(all_sort=all((!is.unsorted(start, strictly = T)) & (!is.unsorted(stop, strictly = T)))) %>%
              dplyr::summarise(as=all(all_sort)) %>%
              dplyr::pull(as),
            all(break_df$start<break_df$stop)
  )
  return(dplyr::mutate(snp_df, region_id=set_ld_region(break_df, snp_df, assign_all = assign_all)))
}

assign_region.NULL <- function(break_coord, snp_coord, assign_all){
  liftover_allf <- system.file("fourier_ls-all.bed.gz", package="ldshrink")
  break_df <- readr::read_delim(liftover_allf, delim="\t", trim_ws = T, skip=1, col_names=c("chr", "start", "stop"))
  break_df <- dplyr::mutate(break_df, chr=as.integer(gsub("chr", "", chr))) %>% dplyr::mutate(region_id=1:n())
  stopifnot(!is.null(break_df$chr),
            !is.null(break_df$start),
            !is.null(break_df$stop),
            !is.null(break_df$region_id),
            is.integer(break_df$chr),
            is.integer(break_df$start),
            is.integer(break_df$stop),
            dplyr::group_by(break_df, chr) %>% dplyr::summarise(all_sort=all((!is.unsorted(start, strictly = T)) & (!is.unsorted(stop, strictly = T)))) %>% dplyr::summarise(as=all(all_sort)) %>% dplyr::pull(as), 
            all(break_df$start<break_df$stop)
  )
  return(dplyr::mutate(snp_df, region_id = set_ld_region(break_df, snp_df, assign_all = assign_all)))
}

assign_region.default <- function(break_coord, snp_coord, assign_all){
  stop("Don't know how to assign region for break_df of type", class(break_coord)[[1]], call. = FALSE)
}



#' Assign SNPs to LD blocks
#'
#' @param snp_df dataframe of snp coordinates, must contain columns named `chr` and `pos`.
#' @param break_df dataframe of LD blocks, if `NULL`, use precomputed (EUR) ldshrink LD blocks, must contain columns named `chr` `start`, `stop`, and `region_id`
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#'
#' @return modified `snp_df` dataframe with additional column `region_id` mapping snp to LD block
#' @export
chunk_genome <- function(snp_df, n_chunks=NA, chunk_size=NA, min_size=10){
    stopifnot(!all(is.na(c(n_chunks, chunk_size))), !all(!is.na(c(n_chunks, chunk_size))))

    if(!is.na(n_chunks)){
        snp_df <- dplyr::group_by(snp_df, chr) %>%
            dplyr::mutate(t_region_id=as.integer(gl(n = n_chunks, k = ceiling(n()/n_chunks), length=n()))) %>% dplyr::ungroup()
    }else{
        snp_df <- dplyr::group_by(snp_df, chr) %>%
            dplyr::mutate(t_region_id=as.integer(gl(n = ceiling(n()/chunk_size) , k = chunk_size, length=n()))) %>% dplyr::ungroup()
        snp_df <- dplyr::group_by(snp_df, chr, t_region_id) %>% dplyr::summarise(ct=n()) %>% dplyr::ungroup() %>% dplyr::inner_join(snp_df) %>% dplyr::mutate(t_region_id=ifelse(ct<min_size, t_region_id-1, t_region_id)) %>% dplyr::select(-ct) %>% dplyr::ungroup()
    }
    snp_df <- dplyr::distinct(snp_df, chr, t_region_id) %>% dplyr::mutate(region_id=1:n()) %>% dplyr::inner_join(snp_df) %>% dplyr::select(-t_region_id)
    return(snp_df)
}



#' Assign or Interpolate Genetic Map Values
#'
#' @param snp_df 
#' @param map_df dataframe of reference genetic map values
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#' @return a modified `snp_df` with an additional column giving the interpolated genetic map values at each locus
#' @export
assign_genetic_map <- function(snp_df, map_df, strict=FALSE){
  u_chr <- dplyr::distinct(snp_df, chr)
  snp_dfl <- split(snp_df, snp_df$chr)
  map_dfl <- dplyr::semi_join(map_df, u_chr, by="chr") %>% split(.$chr)
  stopifnot(all(names(map_dfl)==names(snp_dfl)))
  retdf <- purrr::map2_df(map_dfl, snp_dfl, ~dplyr::mutate(.y, map=interpolate_genetic_map(.x$map, .x$pos, .y$pos)))
  return(retdf)
}


calc_theta <- function(m){
  nmsum <- sum(1 / (1:(2*m-1)))
  (1/nmsum) / (2*m + 1/nmsum)
}

flip_allele_exp <- function(allele_a, allele_b){
  utf_i <- Vectorize(utf8ToInt)
  data_df <- tibble::data_frame(allele_a = allele_a, allele_b=allele_b)
  gwas_snp_df <- data_df %>%
    tidyr::separate(allele_a, c("ref_a", "alt_a")) %>%
    tidyr::separate(allele_b, c("ref_b", "alt_b")) %>%
    dplyr::mutate(
      ref_a=utf_i(tolower(ref_a)),
      alt_a=utf_i(tolower(alt_a)),
      ref_b=utf_i(tolower(ref_b)),
      alt_b=utf_i(tolower(alt_b))
    )
 
  return(flip_allele(gwas_snp_df$ref_a,
                     gwas_snp_df$alt_a,
                     gwas_snp_df$ref_b,
                     gwas_snp_df$alt_b))
}
