library(tidyverse)



map_files <- snakemake@input[["mapf"]]
outf <- snakemake@output[["mapf"]]


map_file_df <- data_frame(map_file = map_files,
                          chr = gsub(".*chr([0-9]+).+",
                                     "\\1", basename(map_files)))

map_df <- rowwise(map_file_df) %>% do(map2_dfr(.$map_file, .$chr,
                                               function(filename, chrom){
                                                   read_delim(filename,
                                                              col_names = c("SNP", "pos", "map"),
                                                              delim = " ") %>% mutate(chr = chrom)
                                               }))


saveRDS(map_df, outf)
