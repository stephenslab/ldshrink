library(SeqArray)
library(tidyverse)
library(LDshrink)



outf <- snakemake@output[["gdsf"]]
toutf <- snakemake@input[["temp_gds"]]
breakf <- snakemake@input[["breakf"]]



save.image()
import_panel_data(temp_gds=toutf,
                  output_file = outf,
                  ld_break_file=breakf,
                  overwrite=F,
                  parallel=cores)
                                                

