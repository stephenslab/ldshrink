library(LDshrink)
library(tidyverse)
library(SeqArray)
library(RSSp)


evdf <- snakemake@input[["evd_chunkf"]]
gdsf <- snakemake@input[["gdsf"]]
chrom <- as.character(snakemake@params[["chrom"]])
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))
outf <- snakemake@output[["outf"]]
soutf <- snakemake@output[["soutf"]]


save.image()

Ql_df <- read_si_ql(evdf)
Ql <- Ql_df[["Ql"]]
Dl <- Ql_df[["Dl"]]



gds <- seqOpen(gdsf, readonly = T)
seqSetFilterChrom(gds, chrom)

n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))

res_l <- gen_sim_gds_direct_ldsc(Ql = Ql, Dl = Dl, gds = gds,
                                 pve = pve, bias = bias,
                                 nreps = nreps, fgeneid = mfgeneid)

stopifnot(all(sort(names(res_l$ldsc_df_l)) == sort(mfgeneid)))
tparam_df <- mutate(gen_tparamdf_norm(pve = pve,bias = bias,nreps = nreps,n = n,p = p,fgeneid = mfgeneid),
                    n = n, p = p)

for (i in names(res_l$ldsc_df_l)){
    stopifnot(nrow(res_l$ldsc_df_l[[i]]) == p)
    write_delim(res_l$ldsc_df_l[[i]], path = outf[as.integer(i)], delim = "\t")
    ttparam_df <- tparam_df %>% filter(fgeneid == i)
    write_delim(ttparam_df, path = soutf[as.integer(i)], delim = "\t")
}
