library(LDshrink)
library(RcppEigenH5)
library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))
cores <- as.integer(as.numeric(snakemake@threads))
LDchunk <- snakemake@params[["LDchunk"]]
LDchunkf <- snakemake@output[["LDchunkf"]]

rdsf <- snakemake@output[["rdsf"]]



gds <- seqOpen(gdsf, readonly = T)
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
                     nreps = nreps, fgeneid = c(mfgeneid),cores=cores)

LD_region <- seqGetData(gds, "annotation/info/LD_chunk")
p <- length(LD_region)
split_LDR <- split(1:p,LD_region)
stopifnot(length(split_LDR)==length(LDchunk))
pb <- progress_bar$new(total=length(split_LDR))
for(i in 1:length(split_LDR)){
    tldchunk <- names(split_LDR)[i]
    j <- which(tldchunk==LDchunk)
    write_mat_h5(h5file = LDchunkf[j],groupname = tldchunk,dataname = "uh",data = res_l$bias_uh_mat[split_LDR[[i]],],deflate_level = 0,doTranspose = F)
    pb$tick()
}
saveRDS(res_l,rdsf)

