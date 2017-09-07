library(SeqArray)
library(tidyverse)
library(LDshrink)


inpf <- snakemake@input[["vcff"]]
outf <- snakemake@output[["gdsf"]]
toutf <- snakemake@output[["temp_gds"]]
mapf <- snakemake@input[["mapf"]]
panelf <- snakemake@input[["panelf"]]
breakf <- snakemake@input[["breakf"]]

cores <- as.integer(snakemake@threads)


map_df <- readRDS(mapf)
import_panel_data(vcf_files = inpf,
                  map_df=map_df,
                  output_file = outf,
                  temp_gds=toutf,
                  overwrite=F,
                  parallel=cores,
                  panel_ind_file=panelf,
                  subset_MAF=0.05)
                                                

