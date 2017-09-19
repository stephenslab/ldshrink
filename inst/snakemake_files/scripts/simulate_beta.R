library(LDshrink)
library(RcppEigenH5)
library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
snp_df <- readRDS(snakemake@input[["snp_infof"]])
pve <- as.numeric(snakemake@params[["pve"]])

out_rdsf <- snakemake@output[["beta_dff"]]

gds <- seqOpen(gdsf, readonly = T)
subset_gds(gds,snp_df)

n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))




tparam_df <- data_frame(
    tpve = pve,
    tsigu = sqrt(n / p * tpve)) %>%
    mutate(n = n, p = p)

beta_df <- sim_beta_gds(
    gds,
    tparam_df = tparam_df
)

saveRDS(beta_df,out_rdsf)
