library(RcppEigenH5)


evdf <- snakemake@input[["evdf"]]
uhf <- snakemake@input[["uhf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
quhf <- snakemake@input[["quhf"]]


Q_chunk <- read_2d_mat_h5(evdf,LDchunk,"Q")
uh_chunk <- read_2d_mat_h5(uhf,LDchunk,"uh")

write_mat_h5(quhf,LDchunk,"quh",crossprod(Q_chunk,uh_chunk),deflate_level=0,doTranspose=F)

