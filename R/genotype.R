

parse_listop <- function(x, names){
    stopifnot(length(x)>=1, length(names)>1)
    purrr::reduce(names, function(x, y, z){ z[[x]] %||% z[[y]] }, z=x)
}





#' R6 object for representing refernce panel data
#'
#' the `reference_panel` class encapsulates all the information necessary for computing LD from a reference panel
#'
## reference_panel <- R6Class("reference_panel",
##                            public = list(
##                                data = NULL,
##                                snp_id = NULL,
##                                sample_id = NULL,
##                                is_haplotype = FALSE,
##                                chrom=NULL,
##                                pos=NULL,
##                                map=NULL,
##                                ref=NULL,
##                                alt=NULL,
##                                population=NULL,
##                                snp_index=NULL,
##                                sample_index=NULL,
##                                options=NULL,
##                                initialize = function(panel_matrix, snpinfo=list(), ...){
##                                    self$data  <-  panel_matrix
##                                    argl <- list(...)
##                                    p <- nrow(snpinfo)
##                                    stopifnot(ncol(self$data)==p)
##                                    self$is_haplotype <- argl[["is_haplotype"]] %||% dplyr::between(range(c(self$data)), 0, 1)
##                                    snp_id_opts <- c("snp_id", "rsid", "SNP", "snp","variant_id")
##                                    chrom_opts <-c("chrom", "chr", "chromosome")
##                                    allele_opts <-c("allele", "ALLELE")
##                                    ref_opts <- c("ref", "a1", "A1", "REF","non_effect_allele")
##                                    alt_opts <- c("alt", "a2", "A2", "ALT","effect_allele")
##                                    pos_opts <- c("pos", "position", "POS")
##                                    map_opts <- c("map", "genetic_map")
##                                    all_opts <- c(
##                                        snp_id_opts,
##                                        chrom_opts,
##                                        ref_opts,
##                                        alt_opts,
##                                        pos_opts,
##                                        map_opts, "snp_index", "sample_index")
##                                    self$snp_id <- parse_listop(snpinfo, snp_id_opts) %||%
##                                        parse_listop(argl, snp_id_opts) %||% colnames(panel_matrix) %||%
##                                        1L:p
##                                    self$snp_index <- argl[["snp_index"]] %||% self$snp_id

##                                    self$chrom <- parse_listop(snpinfo, chrom_opts) %||% parse_listop(argl, chrom_opts)



##                                    ref <- parse_listop(snpinfo, ref_opts) %||% parse_listop(argl, ref_opts)
##                                    alt <- parse_listop(snpinfo, alt_opts) %||% parse_listop(argl, alt_opts)
##                                    if(is.null(ref) ||is.null(alt)){
##                                        stop("No allele information provided, `allele` stored as `NULL`")
##                                    }
##                                    self$pos  <- parse_listop(snpinfo, pos_opts) %||% parse_listop(argl, pos_opts)
##                                    self$map <- parse_listop(snpinfo, map_opts) %||% parse_listop(argl, map_opts)
##                                    self$sample_id <- argl[["sample_id"]] %||% rownames(panel_matrix)
##                                    self$sample_index <- argl[["sample_index"]] %||% self$sample_id
##                                    self$options <- argl[! names(argl) %in% all_opts]
##                                    if(length(names(self$options))>0){
##                                        message(paste0("Keeping options: ", paste0(names(self$options), collapse=" , ")))
##                                    }
##                                }
##                            )
##                            )
