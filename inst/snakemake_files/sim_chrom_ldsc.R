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
outf <- snakemake@output[["outf"]]
soutf <- snakemake@output[["soutf"]]



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



res_l <- gen_sim_gds_ldsc(gds, pve = pve, bias = bias,
                          nreps = nreps, fgeneid = c(mfgeneid))
stopifnot(all(sort(names(res_l$ldsc_df_l))==sort(mfgeneid)))

for (i in names(res_l$ldsc_df_l)){
    stopifnot(nrow(res_l$ldsc_df_l[[i]]) == p)
    write_delim(res_l$ldsc_df_l[[i]], path = outf[as.integer(i)], delim = "\t")
    tparam_df <- mutate(res_l[["tparam_df"]],
                        n = res_l[["n"]], p = res_l[["p"]]) %>%
        filter(fgeneid == i)
    write_delim(tparam_df, path = soutf[as.integer(i)], delim = "\t")
}
