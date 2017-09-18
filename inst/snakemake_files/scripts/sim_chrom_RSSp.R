library(LDshrink)
library(tidyverse)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
chrom <- as.character(snakemake@params[["chrom"]])
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))

rdsf <- snakemake@output[["rdsf"]]



gds <- seqOpen(gdsf, readonly = T)
seqSetFilterChrom(gds, chrom)
n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))
cat("N: ", n, "\n")
cat("PVE: ", pve, "\n")
cat("bias: ", bias, "\n")
tparam_df <- gen_tparamdf_norm(pve, bias, nreps, n, p) %>%
    filter(fgeneid %in% mfgeneid)
mpve <- tparam_df$tpve
mbias <- tparam_df$tbias


save.image()

res_l <- gen_sim_gds(gds, pve = pve, bias = bias,
                          nreps = nreps, fgeneid = c(mfgeneid))

saveRDS(res_l,rdsf)

