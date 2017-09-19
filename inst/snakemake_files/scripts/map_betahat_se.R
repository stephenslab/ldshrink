library(LDshrink)
library(RcppEigenH5)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
outf <- snakemake@output[["sumstatf"]]
ymatf <- snakemake@input[["ymatf"]]
ymat <- read_2d_mat_h5(snakemake@input[["ymatf"]],
                       "trait",
                       "ymat")



gds <- seqOpen(gdsf, readonly = T)
h5grps <- h5ls(ymatf)
if(sum(grepl("SNPinfo",h5grps))){
    snpi_df <- read_df_h5(ymatf,"SNPinfo")
    subset_gds(gds,snpi_df)
}


map_bh_se_gds(gds,ymat,out_file=outf)
