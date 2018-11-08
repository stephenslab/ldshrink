#' Calculate LD from a reference panel
#'
#' `ldshrink` constructs an LD matrix from a reference panel.
#'
#'
#'
#' @param reference_panel is an `n` by `p` matrix of reference-panel haplotypes or genotypes which will be used for LD calculation
#' @param map is a length `p` vector of cumulative genetic map values. `map` must be _strictly_ _sorted_
#' @param method_ld method for computing LD,(defaults to the Wen-Stephens 2010 shrinkage-estimator)  a list of available methods may be found in the `Details` section below
#' @param method_ld_params a list of additional arguments to be passed to the LD method
#' @param output a length one character vector specifying the format for output.  defaults to a sparse,symmetric matrix.
#' @param store_anno a boolean indicating whether to store the result of the distance function, can only be used with `output="data.frame"`
#' @param progress a boolean indicating whether to show a progress bar
#' @param progress isGenotype a boolean indicating whether the data is genotype data (the alternative being that the data is haplotype data)
#' @return panel_ld an estiamate of LD
#' @export
#'
estimate_LD <- function(reference_panel,
                        map,
                        method="ldshrink",
                        method_ld_params=list(m=85, Ne=11490.672741,cutoff=1e-3),
                        output="dsCMatrix",
                        store_anno = FALSE,
                        progress=FALSE,
                        isGenotype=TRUE,
                        ...){

    colnames(reference_panel) <- colnames(reference_panel) %||% (names(map) %||% as.character(1:ncol(reference_panel)))
    if(output=="matrix"){
        toutput <- "dsCMatrix"
    }else{
        toutput <- output
    }
    if(store_anno){
        stopifnot(output %in% c("data.frame","data_frame"))
    }
  argl <- list(...)
    optl <- c(
        list(is_haplotype=!isGenotype,
             write_annotations=store_anno,
             progress=progress,
             output=toutput),
        method_ld_params,
        argl)
    if(method=="ldshrink"){
        ret <- ldshrink_cor(reference_panel,map,optl)
    }else{
        ret <-sample_cor(reference_panel,map,optl)
    }
    if(output=="matrix"){
        return(as.matrix(ret))
    }
    return(ret)
}



ldshrink2 <- function(reference_panel.x, reference_panel.y, method_ld="ldshrink", method_dist="ldshrink",
                      method_ld_params=list(m=85, Ne=11490.672741, isGenotype=TRUE),
                      methold_dist_params=list(m=85, Ne=11490.672741),
                      dist_cutoff = 1e-3,
                      store_anno = FALSE, ...){





}



#' @param reference_panel is an R6 object of type "reference_panel".  (For more info, see `?reference_panel`
#' @param method_ld method for computing LD, a list of available methods may be found in the `Details` section below
#' @param method_dist method for computing the distance betwteen two variants (usually based on physical distance or genetic map distance)
#' @param method_ld_params a list of additional arguments to be passed to the LD method
#' @param method_dist_params a list of additional arguments to be passed to the distance method
#' @param dist_cutoff threshold for the output of `method_dist` below which LD will not be calculated, set to `NULL` will result in LD being calculated for every pair of SNPs
#' @param store_anno a boolean indicating whether to store the result of the distance function
#' @return panel_ld
#' @export
#'
ldshrink_evd <- function(panel, map=NULL, m=85,
                            Ne=11490.672741,
                            cutoff=1e-3,
                         useldshrink=T, na.rm=F,...){
    isGeno <- max(panel,na.rm = na.rm)>1
    if(useldshrink){
        panel  <-  scale(panel,center=TRUE,scale=FALSE)
        S <-estimate_LD(reference_panel = panel,
                        map = map,
                        method="ldshrink",
                        method_ld_params=list(m = m,Ne = Ne,cutoff = cutoff),
                        output="matrix",...)
    }else{
      S <- stats::cor(panel, use = "complete.obs")
    }
    N <- nrow(panel)
    L2 <- estimate_LDscores(S)
    evdR <- eigen(S)
    return(list(R=S, L2=L2, D=evdR$values, Q=evdR$vectors))
}



#' Estimate LD scores from the LD matrix
#'
#' @param R a pxp LD matrix (as obtained from `estimate_LD`)
#' @param N a scalar representing the number of samples in the reference panel
#'
#' @return A length p vector with LD scores at each locus
#' @export
#'
#' @examples
estimate_LDscores <- function(R,N){
  L2 <- colSums(S^2)-1
  L2-(1-L2)/(N-2)
}
