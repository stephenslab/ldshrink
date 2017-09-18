library(LDshrink)
library(RcppEigenH5)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
outf <- snakemake@output[["sumstatf"]]
ymat <- read_2d_mat_h5(snakemake@input[["ymatf"]],
                       "trait",
                       "ymat")


gds <- seqOpen(gdsf, readonly = T)


map_bh_se_gds(gds,ymat,out_file=outf)
