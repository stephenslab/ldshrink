library(SeqArray)

gdsf <- snakemake@input[["gdsf"]]


gds <- seqOpen(gdsf,readonly = T)
LD_chunks <- unique(seqGetData(gds,var.name="annotation/info/LD_chunk"))
for(i in LD_chunks){
    write(i,file=paste0(i,".txt"))
}
seqClose(gds)

