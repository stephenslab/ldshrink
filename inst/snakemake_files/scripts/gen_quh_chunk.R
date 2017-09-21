library(RcppEigenH5)
library(RSSp)

evdf <- snakemake@input[["evdf"]]
uhf <- snakemake@input[["uhf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
quhf <- snakemake@output[["quhf"]]



gw_snpi <- read_ivec(uhf,LDchunk,"snp_id")
uh_chunk <- read_2d_mat_h5(uhf, LDchunk,"uh")
stopifnot(!any(is.na(uh_chunk)))
resl <- gen_quh_chunk_mat(uh_chunk,evdf,gw_snpi)

write_mat_h5(quhf,
             LDchunk,             
             "quh",
             resl$quh,
             deflate_level=0,
             doTranspose=F)

write_dvec_h5(quhf,LDchunk,"D",resl$D)

write_ivec_h5(quhf,LDchunk,"D",resl$snp_id)
