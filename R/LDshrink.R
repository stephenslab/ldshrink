
approx_diag_l <- function(R, indl) {
    return(lapply(indl, function(x, R) {
        return(R[x, x])
    }, R = R))
}


make_chunks <- function(R, chunksize) {
    library(BBmisc)
    p <- nrow(R)
    stopifnot(p == ncol(R))
    chunkl <- chunk(1:p, chunk.size = chunksize)
    diag_Rl <- approx_diag_l(R, chunkl)
    return(list(index = chunkl, matrices = diag_Rl))
}


el_or_default <- function(el, defaults, name){
  if (is.null(el[[name]])) {
    return(defaults[[name]])
  }else{
    return(el[[name]])
  }
}


make_chunks <- function(R, chunksize = ncol(R)) {
    library(BBmisc)
    p <- nrow(R)
    stopifnot(p == ncol(R))
    chunkl <- chunk(1:p, chunk.size = chunksize)
    diag_Rl <- approx_diag_l(R, chunkl)
    return(list(index = chunkl, matrices = diag_Rl))
}

write_eigen_chunks <- function(genof = "", evdf = "", chunksize, dataname = "R", R = NULL, write = T, map = NULL) {
    library(progress)
    library(RcppEigenH5)
    if (file.exists(evdf)) {
        if (write) {
            stopifnot(!group_exists(evdf, as.character(chunksize)))
        }
    }
    stopifnot(!is.null(R) | file.exists(genof))
    if (is.null(R)) {
        R <- read_2d_mat_h5(genof, groupname = "/", dataname = dataname)
    }
    mchunks <- make_chunks(R, chunksize)
    evdRl <- lapply(mchunks[["matrices"]], eigen, symmetric = T)
    Ql <- list()
    rdl <- list()
    indl <- list()
    i <- 1
    nchunks <- length(evdRl)
    pb <- progress_bar$new(total = nchunks)
    rm(R)
    gc()
    for (i in 1:nchunks) {
        Ql[[i]] <- evdRl[[i]]$vectors
        rdl[[i]] <- evdRl[[i]]$values
        indl[[i]] <- mchunks[["index"]][[i]]
        if (write) {
            write_dvec_h5(evdf, paste0(chunksize, "/", i), "D", rdl[[i]], deflate_level = 4)
            write_mat_h5(evdf, paste0(chunksize, "/", i), "Q", data = Ql[[i]], deflate_level = 2)
            write_ivec_h5(evdf, paste0(chunksize, "/", i), "ind", data = indl[[i]], deflate_level = 4)
        }
        pb$tick()
    }
    return(list(Ql = Ql, rdl = rdl, indl = indl))
}





read_eigen_chunks <- function(evdf, chunksize) {
    # library(RcppEigenH5)
    stopifnot(file.exists(evdf))
    stopifnot(RcppEigenH5::group_exists(evdf, as.character(chunksize)))
    all_chunks <- RcppEigenH5::list_groups_h5(evdf, as.character(chunksize))
    chunk_grps <- paste0(chunksize, "/", all_chunks)
    Ql <- list()
    rdl <- list()
    indl <- list()
    for (i in 1:length(chunk_grps)) {
        chunk_grp <- chunk_grps[i]
        Ql[[i]] <- RcppEigenH5::read_2d_mat_h5(h5file = evdf, groupname = chunk_grp, dataname = "Q")
        rdl[[i]] <- RcppEigenH5::read_dvec(h5file = evdf, groupname = chunk_grp, dataname = "D")
        indl[[i]] <- RcppEigenH5::read_ivec(h5file = evdf, groupname = chunk_grp, dataname = "ind")
    }
    names(Ql) <- chunk_grps
    names(rdl) <- chunk_grps
    names(indl) <- chunk_grps
    return(list(Ql = Ql, Dl = rdl, indl = indl))
}


#' Chunk LD as in ldetect
#' Given a genetic map dataframe, creates LD chunks
#' @param map_dat dataframe with a column 'map' 
#' with cumulative genetic map data, and a column 'pos'
#'  with the physical coordinate of each SNP
chunk_LD <- function(map_dat, panel_size = 379, N = 5000,
                     Ne = 11418, cutoff = 1.5e-08, progress = F) {
    
  library(progress)

    nsnp <- nrow(map_dat)
    chunknum <- floor(nsnp/N)
    nind <- panel_size
    if (progress) {
      pb <- progress_bar$new(total = chunknum)
    }    
    map <- map_dat$map
    
    pos <- map_dat$pos
    ipos <- 1:nsnp
    
    startpos <- integer(chunknum)
    startipos <- integer(chunknum)
    endpos <- integer(chunknum)
    endipos <- integer(chunknum)
    for (i in 1:chunknum) {
        start <- (i - 1) * N + 1
        end <- (i - 1) * N + N + 1
        if (i == chunknum) {
            end <- length(pos)
            startpos[i] <- pos[start]
            startipos[i] <- start
            endpos[i] <- pos[end]
            endipos[i] <- end
            break
        }
        startpos[i] <- pos[start]
        startipos[i] <- start
        endpos[i] <- pos[end - 1]
        endipos[i] <- pos[end - 1]
        endgpos <- map[end - 1]
        test <- end + 1
        testpos <- pos[test]
        stop <- F
        while (!stop) {
            if (test == length(pos)) {
                stop <- T
                # cat('break!\n')
                break
            }
            testpos <- pos[test]
            testgpos <- map[test]
            df <- testgpos - endgpos
            tmp <- exp(-4 * Ne * df/(2 * nind))
            if (tmp < cutoff) {
                stop <- T
            } else {
                test <- test + 1
            }
        }
        endipos[i] <- test
        endpos[i] <- testpos
        if (progress) {
          pb$tick()
        }
    }
    return(tibble::data_frame(startpos = startpos, endpos = endpos,
                              startipos = startipos, endipos = endipos))
}

