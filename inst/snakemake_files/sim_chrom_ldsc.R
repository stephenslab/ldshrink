library(LDshrink)
library(tidyverse)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
chrom <- as.character(snakemake@params[["chrom"]])
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
res_l <- readRDS(snakemake@input[["rdsf"]])
outf <- snakemake@output[["outf"]]
soutf <- snakemake@output[["soutf"]]
rdsf <- snakemake@output[["rdsf"]]



gds <- seqOpen(gdsf, readonly = T)
seqSetFilterChrom(gds, chrom)
# save.image()
# stop()
res_l[["ldsc_df_l"]] <- read_SNPinfo_ldsc_gwas(gds,res_l$bias_uh_mat,N=res_l$n)
n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))

# res_l <- gen_sim_gds_ldsc(gds, pve = pve, bias = bias,
#                           nreps = nreps, fgeneid = c(mfgeneid))

# saveRDS(res_l,rdsf)
stopifnot(all(sort(names(res_l$ldsc_df_l))==sort(mfgeneid)))

# 
# all_rssp <- prep_RSSp_evd(Ql=ql_df$Ql,Dl=ql_df$Dl,U = res_l$bias_uh_mat,N = res_l$n)
# est_sim()

for (i in names(res_l$ldsc_df_l)){
    stopifnot(nrow(res_l$ldsc_df_l[[i]]) == p)
    write_delim(res_l$ldsc_df_l[[i]], path = outf[as.integer(i)], delim = "\t")
    tparam_df <- mutate(res_l[["tparam_df"]],
                        n = res_l[["n"]], p = res_l[["p"]]) %>%
        filter(fgeneid == i)
    write_delim(tparam_df, path = soutf[as.integer(i)], delim = "\t")
}
