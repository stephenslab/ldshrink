library(LDshrink)
library(tidyverse)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
osnpf <- snakemake@output[["snpif"]]
gds <- seqOpen(gdsf)
snpi <- read_SNPinfo_gds(gds)
saveRDS(snpi,osnpf)
