library(RcppEigenH5)


evdf <- snakemake@input[["evdf"]]
uhf <- snakemake@input[["uhf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
quhf <- snakemake@output[["quhf"]]

                                        #save.image()
D_chunk <- read_dvec(evdf,"EVD","D")
Q_chunk <- read_2d_mat_h5(evdf, "EVD","Q")
uh_chunk <- read_2d_mat_h5(uhf, LDchunk,"uh")
prod <- crossprod(Q_chunk, uh_chunk)
write_mat_h5(quhf,             LDchunk,
             "quh",
             prod,
             deflate_level=0,
             doTranspose=F)

i_chunk <- read_ivec(uhf, LDchunk, "snp_id")
write_ivec_h5(quhf, LDchunk, "snp_id",i_chunk)

write_dvec_h5(quhf,LDchunk,"D",D_chunk)
