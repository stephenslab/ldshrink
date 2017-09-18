library(dplyr)
library(RcppEigenH5)
library(purrr)
library(progress)

inf <- snakemake@input[["evdf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
paramf <- as.character(snakemake@input[["phenof"]])
outf <- snakemake@output[["rdsf"]]

nchunks <- length(inf)
iresl <- list()
iresl[["tparam_df"]] <- read_df_h5(paramf,"SimulationInfo")
iresl[["quh_mat"]] <- concat_mat_chunks(inf, LDchunk,
                                        rep("quh", length(LDchunk)))
iresl[["D"]] <- unlist(map2(inf, LDchunk, read_dvec, dataname = "D"))
    
saveRDS(iresl, outf)
