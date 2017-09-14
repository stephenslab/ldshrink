library(dplyr)
library(RcppEigenH5)
library(purrr)
library(progress)

inf <- snakemake@input[["evdf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
iresl <- readRDS(snakemake@input[["rdsf"]])
outf <- snakemake@output[["rdsf"]]

nchunks <- length(inf)
iresl[["quh_mat"]] <- concat_mat_chunks(inf, LDchunk,
                                        rep("quh", length(LDchunk)))
iresl[["D"]] <- unlist(map2(inf, LDchunk, read_dvec, dataname = "D"))
    
saveRDS(iresl, outf)
