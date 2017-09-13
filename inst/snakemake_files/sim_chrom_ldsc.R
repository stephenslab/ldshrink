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


save.image()
# stop()
res_l <- gen_sim_gds_ldsc(gds, pve = pve, bias = bias,
                          nreps = nreps, fgeneid = c(mfgeneid))
stopifnot(all(sort(names(res_l$ldsc_df_l))==sort(mfgeneid)))


#Generate a large block-diagonal matrix
read_si_ql <- function(chunk_evdf){
  library(RcppEigenH5)
  library(tidyr)
  library(dplyr)

  retl <- list()
  retl[["LDinfo"]] <- read_df_h5(chunk_evdf,"LDinfo") %>% nest(-region_id)
  region_id <-  retl[["LDinfo"]][["region_id"]]
  retl[["Dl"]] <- lapply(region_id,read_dvec,h5file=chunk_evdf,dataname="D")
  
  names(retl[["Dl"]]) <- region_id
  retl[["Ql"]] <- lapply(region_id,read_2d_mat_h5,h5file=chunk_evdf,dataname="Q")
  names(retl[["Ql"]]) <- region_id
  return(retl)
}
all_rssp <- prep_RSSp_evd(Ql=ql_df$Ql,Dl=ql_df$Dl,U = res_l$bias_uh_mat,N = res_l$n)
est_sim()

for (i in names(res_l$ldsc_df_l)){
    stopifnot(nrow(res_l$ldsc_df_l[[i]]) == p)
    write_delim(res_l$ldsc_df_l[[i]], path = outf[as.integer(i)], delim = "\t")
    tparam_df <- mutate(res_l[["tparam_df"]],
                        n = res_l[["n"]], p = res_l[["p"]]) %>%
        filter(fgeneid == i)
    write_delim(tparam_df, path = soutf[as.integer(i)], delim = "\t")
}
