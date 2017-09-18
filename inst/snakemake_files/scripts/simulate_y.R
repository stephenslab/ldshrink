library(LDshrink)
library(RcppEigenH5)
library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
betamat <- readRDS(snakemake@input[["betamatf"]])
pve <- as.numeric(snakemake@params[["pve"]])
tgdsf <- tempfile()
file.copy(gdsf,tgdsf)
if(ncol(betamat)!=length(pve))
    stopifnot(length(pve)==1)

out_h5f <- snakemake@output[["h5f"]]
gds <- seqOpen(tgdsf, readonly = F)

n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))
stopifnot(p == nrow(betamat))
add.gdsn(index.gdsn(gds, "/annotation/info/"), "beta_mat", betamat)
tparam_df <- data_frame(tpve = pve, tsigu = sqrt(n / p * tpve)) %>%
    mutate(n = n, p = p)

ymat <- gen_sim_phenotype_gds(gds, tparam_df = tparam_df, fgeneid = mfgeneid,cores = cores)
write_mat_h5(out_h5f, "trait", "ymat", ymat, deflate_level = 0L)
write_df_h5(tparam_df, groupname = "SimulationInfo", outfile = out_h5f, deflate_level = 0L)
