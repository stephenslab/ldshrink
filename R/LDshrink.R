#' Calculate LD from a reference panel
#'
#' `ldshrink` constructs an LD matrix from a reference panel.
#'
#'
#'
#' @param reference_panel is an `n` by `p` matrix of reference-panel haplotypes or genotypes which will be used for LD calculation
#' @param map is a length `p` vector of cumulative genetic map values. `map` must be _strictly_ _sorted_
#' @param method method for computing LD,(defaults to the Wen-Stephens 2010 shrinkage-estimator)  a list of available methods may be found in the `Details` section below
#' @param method_ld_params a list of additional arguments to be passed to the LD method
#' @param output a length one character vector specifying the format for output.  defaults to a sparse,symmetric matrix.
#' @param store_anno a boolean indicating whether to store the result of the distance function, can only be used with `output="data.frame"`
#' @param progress a boolean indicating whether to show a progress bar
#' @param isGenotype a boolean indicating whether the data is genotype data (the alternative being that the data is haplotype data)
#' @return panel_ld an estiamate of LD
#' @export
#' @examples
#' #load data
#' data(reference_genotype)
#' data(reference_map)
#'#obtain a symmetric sparse
#' sparse_LD  <- estimate_LD(reference_genotype,reference_map)
#' #obtain a dense matrix
#' dense_LD  <- estimate_LD(reference_genotype,reference_map,output="matrix")
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



#' Calculate LD and eigenvalue decomposition from a reference panel
#'
#' @param reference_panel reference panel haplotype
#' @param map refrence genetic map
#' @param m sample size for genetic map
#' @param Ne effective population size
#' @param cutoff correlation cutoff
#' @param useldshrink whether or not to use LD shrink
#' @param na.rm whether or not to remove missing values
#' @param ... extra arguments to other functions(
#' @export
ldshrink_evd <- function(reference_panel, map=NULL, m=85,
                            Ne=11490.672741,
                            cutoff=1e-3,
                         useldshrink=T, na.rm=F, ...){
    isGeno <- max(reference_panel,na.rm = na.rm) > 1
    if (useldshrink) {
        # panel  <-  scale(panel,center=TRUE,scale=FALSE)
        S <- estimate_LD(reference_panel = reference_panel,
                        map = map,
                        method = "ldshrink",
                        method_ld_params = list(m = m, Ne = Ne, cutoff = cutoff),
                        output = "matrix",isGenotype = isGeno,...)
    }else{
      S <- stats::cor(reference_panel, use = "complete.obs")
    }
    N <- nrow(reference_panel)
    L2 <- estimate_LDscores(S,N)
    evdR <- eigen(S)
    return(list(R = S, L2 = L2, D = evdR$values, Q = evdR$vectors))
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
#' data(reference_genotype)
#' data(reference_map)
#' R <- estimate_LD(reference_panel=reference_genotype,map=reference_map)
#' L2 <- estimate_LDscores(R,nrow(reference_genotype))
estimate_LDscores <- function(R,N){
  denom <- ifelse(N>2,N-2,N)
    if(inherits(R,"Matrix")||inherits(R,"matrix")){
        L2 <- apply(R^2-(1-R^2)/denom,2,sum)
    }else{
        if(inherits(R,"data.frame")){
            L2 <- purrr::map_dbl(R$data,~sum(.x$r^2-(1-.x$r^2)/denom))-1
        }else{
            stop("R must be a data.frame,matrix or Matrix")
        }
    }
  return(L2)
}
