
library(SeqArray)
library(LDshrink)
library(readr)
library(dplyr)
library(RSSp)
library(purrr)

inputf <- snakemake@input[["inputf"]]
gdsf <- snakemake@input[["geno_gdsf"]]
ogdsf <- snakemake@output[["geno_gdsf"]]
trait_name <-snakemake@params[["traitn"]]
out_dir <- paste0("../chunks_",trait_name)
outf <- snakemake@output[["output_gwasf"]]


gds <- seqOpen(gdsf)



input_df <- read_delim(inputf,delim="\t")



output_df <- inner_join(input_df,subset_gds(gds,input_df,region_id=T))%>% select(snp_id,Z,N,region_id)

p <- calc_p(gds)

stopifnot(nrow(output_df)==p)

write_delim(output_df,outf,delim="\t",col_names=T)

outl <- split(output_df,output_df$region_id)

iwalk(outl,function(df,name,dirname){
    write_delim(df,file.path(dirname,paste0(name,"_ldsc_sub.txt")),delim="\t")
    },dirname=out_dir)


