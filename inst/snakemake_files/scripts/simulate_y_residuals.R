library(LDshrink)
library(RcppEigenH5)
library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
beta_df <- readRDS(snakemake@input[["beta_dff"]])
residmat <- readRDS(snakemakee@input[["resid_f"]])

pve <- as.numeric(snakemake@params[["pve"]])
tgdsf <- tempfile()
file.copy(gdsf, tgdsf)
if(ncol(betamat) != length(pve))
    stopifnot(length(pve)==1)

out_h5f <- snakemake@output[["h5f"]]
gds <- seqOpen(tgdsf, readonly = F)

subset_gds(gds,beta_df)
index_cols <- colnames(beta_df)
index_cols <- index_cols[! index_cols %in% c("fgeneid","beta")]

betamat <- spread(key = fgeneid,value = beta) %>% select(-one_of(index_cols)) %>% data.matrix
n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))
stopifnot(p == nrow(betamat))


add.gdsn(index.gdsn(gds, "/annotation/info/"), "beta_mat", betamat)
tparam_df <- data_frame(tpve = pve, tsigu = sqrt(n / p * tpve)) %>%
    mutate(n = n, p = p)

ymat <- gen_sim_phenotype_gds(
    gds,
    tparam_df = tparam_df,
    fgeneid = mfgeneid,
    use_beta = T,
    residmat = residmat
)
write_mat_h5(out_h5f, "trait", "ymat", ymat, deflate_level = 0L)
write_df_h5(tparam_df, groupname = "SimulationInfo", outfile = out_h5f, deflate_level = 0L)
