library(SeqArray)

inpf <- snakemake@input[["vcff"]]
outf <- snakemake@output[["gdsf"]]
cores <- as.integer(snakemake@threads)

seqVCF2GDS(inpf, outf, parallel = cores,storage.option="LZ4_RA.fast")

