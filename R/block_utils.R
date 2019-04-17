#' Assign SNPs to LD blocks
#'
#' @param snp_chr dataframe of snp coordinates, must contain (integer-valued) columns named `chr` and `pos`.
#' @param break_chr integer vector of chromosome assignment for each LD block
#' @param break_start integer vector of the start coordinate for each LD block
#' @param break_stop integer vector of the stop coordinate for each LD block
#' @param break_id integer vector with a unique id for each LD block
#' @param snp_chr integer vector giving the chromosome of every SNP
#' @param snp_pos integer vector giving the coordinate of each SNP
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#' @return modified `snp_df` dataframe with additional column `region_id` mapping snp to LD block
#' @export
assign_region <- function(break_chr, break_start, break_stop, break_id = seq_along(break_chr), snp_chr, snp_pos, assign_all=T){
  stopifnot(!is.null(break_chr),
            length(break_chr)==length(break_start),
            length(break_start)==length(break_stop),
            length(break_stop)==length(break_id),
            is.integer(break_start),
            is.integer(break_stop))
  if(is.character(break_chr)|is.character(snp_chr)){
    break_chr <- factor(break_chr,levels = unique(break_chr))
    snp_chr <- factor(snp_chr,levels = levels(break_chr))
  }
  return(set_ld_region(ld_chr = break_chr,
                       ld_start = break_start,
                       ld_stop = break_stop,
                       ld_region_id = break_id,
                       chr = snp_chr,
                       pos = snp_pos,
                       assign_all = assign_all))
}


#' Assign SNPs to LD blocks
#'
#' @param snp_df dataframe of snp coordinates, must contain (integer-valued) columns named `chr` and `pos`.
#' @param n_chunks number of (approximately) equally sized chunks to be returned. (cannot be combined with `chunk_size`)
#' @param chunk_size maximum size of each chunk (cannot be combined with `n_chunks`)
#' @return modified `snp_df` dataframe with additional column `region_id` mapping snp to LD block
#' @export
chunk_genome <- function(snp_df, n_chunks=NA, chunk_size=NA,min_size=1){
  
    stopifnot(!all(is.na(c(n_chunks, chunk_size))), !all(!is.na(c(n_chunks, chunk_size))))

    if(!is.na(n_chunks)){
        snp_df <- dplyr::group_by(snp_df, chr) %>%
            dplyr::mutate(t_region_id=as.integer(gl(n = n_chunks, k = ceiling(n()/n_chunks), length=n()))) %>% dplyr::ungroup()
    }else{
        snp_df <- dplyr::group_by(snp_df, chr) %>%
            dplyr::mutate(t_region_id=as.integer(gl(n = ceiling(n()/chunk_size) ,
                                                    k = chunk_size, length=n()))) %>%
          dplyr::ungroup()
        snp_df <- dplyr::group_by(snp_df, chr, t_region_id) %>%
          dplyr::summarise(ct=n()) %>%
          dplyr::ungroup() %>%
          dplyr::inner_join(snp_df) %>%
          dplyr::mutate(t_region_id=ifelse(ct<min_size, t_region_id-1, t_region_id)) %>%
          dplyr::select(-ct) %>% dplyr::ungroup()
    }
    snp_df <- dplyr::distinct(snp_df, chr, t_region_id) %>%
      dplyr::mutate(region_id=1:n()) %>%
      dplyr::inner_join(snp_df) %>% 
      dplyr::select(-t_region_id)
    return(snp_df)
}



#' Assign or Interpolate Genetic Map Values
#'
#' @param snp_df dataframe of snp coordinates. `snp_df` must contain (integer-valued) columns named `chr` and `pos`.
#' @param map_df dataframe of reference genetic map values. `map_df` must contain integer-valued columns named `chr` and `pos`,
#' as well as a numeric-valued column called `map`, which must be a _strictly_ sorted vector of cumulative genetic map values.
#' @param strict boolean indicating whether
#' @return a modified `snp_df` with an additional column giving the interpolated genetic map values at each locus
#' @export
assign_genetic_map <- function(snp_df, map_df, strict=FALSE){
  u_chr <- dplyr::distinct(snp_df, chr)
  snp_dfl <- split(snp_df, snp_df$chr)
  map_dfl <- dplyr::semi_join(map_df, u_chr, by="chr") %>% split(.$chr)
  stopifnot(all(names(map_dfl)==names(snp_dfl)))
  retdf <- purrr::map2_df(map_dfl, snp_dfl, ~dplyr::mutate(.y, map=interpolate_genetic_map(.x$map, .x$pos, .y$pos,strict=strict)))
  return(retdf)
}


calc_theta <- function(m){
  nmsum <- sum(1 / (1:(2*m-1)))
  (1/nmsum) / (2*m + 1/nmsum)
}


#' Query a a reference panel for a set of SNPs and see if any of them need a sign flip
#'
#' @param query_chr chromosome for query SNPs coded as integer or character (coding must match `panel_chr`)
#' @param query_pos position of query SNPs
#' @param query_allele allele for query SNPs
#' @param query_block optional vector specifiyng a chunking scheme (default is NULL)
#' @param panel_chr chromosome for the reference SNPs coded as integer or character (coding must match `query_chr`)
#' @param panel_pos position of reference SNPs
#' @param panel_allele allele for reference SNPs
#' @param panel_block optional vector specifying a chunking scheme (default is NULL)
#'
#' @return an `integer` vector of length equal to `length(query_chr)` where a `1` indicates that the query position 
#' is in the reference panel, a `0` indicates that it is not, and a `-1` indicates a match
#' , but that a flip is necessary
#' @export
#'
#' @examples
#' 
#' #read in the example data
#' example_bim <- fs::path_package(package = "ldshrink","test_data/reference_genotype.bim")
#' panel_df <- read.table(example_bim,sep="\t",col.names=c("chr","snp","map","pos","ref","alt"))
#' #create "allele" column
#' panel_allele <- paste0(panel_df$ref,",",panel_df$alt)
#' stopifnot(all(nchar(panel_allele)==3))

#' 
query_reference_snpset <- function(query_chr,query_pos,query_allele,query_block=NULL,panel_chr,panel_pos,panel_allele,panel_block=NULL){
  
  stopifnot(is.null(query_block)==is.null(panel_block))
  
  if((is.null(panel_block)&is.null(query_block))){
    join_cols <- c("chr","pos")
  }else{
    join_cols <- c("chr","pos","block")
  }
  if(!is.null(query_block)){
    query_df <- tibble::tibble(chr=query_chr,pos=query_pos,allele=query_allele,block=query_block)
    panel_df <- tibble::tibble(chr=panel_chr,pos=panel_pos,allele=panel_allele,block=panel_block)
   
  }else{
    query_df <- tibble::tibble(chr=query_chr,pos=query_pos,allele=query_allele)
    panel_df <- tibble::tibble(chr=panel_chr,pos=panel_pos,allele=panel_allele)
  }
  ij_df <- dplyr::inner_join(query_df,panel_df,by=join_cols)
  ij_df <- dplyr::mutate(ij_df,allele_match=flip_alleles(allele.x,allele.y))
  ij_df <- dplyr::right_join(ij_df,query_df)
  ij_df <- dplyr::mutate(ij_df,allele_match=dplyr::if_else(is.na(allele_match),0L,allele_match))
  dplyr::pull(ij_df,allele_match)
  
}

