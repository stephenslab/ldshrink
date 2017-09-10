

hap.fn <- snakemake@input[["hapf"]]
leg.fn <- snakemake@input[["legf"]]
sample.fn <- snakemake@input[["filt_f"]]
out.gdsfn <- snakemake@output[["gdsf"]]
chrom= <- snakemake@params[["chrom"]]
compress.geno="LZ4_RA.fast"
compress.annotation="LZ4_RA.fast"

library(readr)
library(dplyr)
library(SeqArray)
library(SNPRelate)
  
  sample.id <- scan(sample.fn,what=character(),sep="\n")
  sample.id <- sapply(sample.id,function(x)c(paste0(x,"-1"),paste0(x,"-2")))
  leg_df <- read_delim(leg.fn,delim=" ") %>% mutate(allele=paste0(a0,",",a1),chrom=chrom)
  geno_mat <- read_delim(hap.fn,delim=" ",col_names = F)
  geno_mat <- data.matrix(geno_mat)
  p <- nrow(geno_mat)
  gc()
  snpgdsCreateGeno(gds.fn = out.gdsfn,genmat = geno_mat,
                   sample.id = sample.id,
                   snp.id = 1:p,
                   snp.rs.id = leg_df$id,
                   snp.chromosome = leg_df$chrom,
                   snp.position = leg_df$position,
                   snp.allele = leg_df$allele,
                   compress.annotation = compress.annotation,
                   compress.geno = compress.geno)
