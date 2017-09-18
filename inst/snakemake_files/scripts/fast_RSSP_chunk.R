library(RcppEigenH5)
library(SeqArray)
library(RSSp)

library(tidyverse)
library(LDshrink)


out_rdsf <- snakemake@output[["logf"]]
input_resl  <- readRDS(snakemake@input[["rdsf"]])
Ql_df <- read_si_ql(snakemake@input[["evd_chunkf"]])




res_rds <- est_sim(input_resl,
                   Ql_df[["Ql"]],
                   D = unlist(Ql_df[["Dl"]]),
                   doConfound = c(T, F),
                   log_params =  F,
                   useGradient =  T)
saveRDS(res_rds, out_rdsf)
