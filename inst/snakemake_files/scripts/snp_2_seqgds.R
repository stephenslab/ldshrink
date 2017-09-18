library(SeqArray)
library(purrr)

walk2(snakemake@input[["temp_gds"]],snakemake@output[["temp_gds"]],function(inf,outf){
    seqSNP2GDS(gds.fn =inf ,out.fn =outf ,storage.option ="LZ4_RA.fast")
    })

cat("Done!")
