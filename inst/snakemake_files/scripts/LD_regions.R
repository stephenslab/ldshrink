

library(SeqArray)


gdsf <- snakemake@input[["gdsf"]]

gds <- seqOpen(gdsf, readonly = T)
LD_region <- unique(seqGetData(gds, "annotation/info/LD_chunk"))
write(LD_region,snakemake@output[["unique_f"]],sep="\n")

