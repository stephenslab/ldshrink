library(readr)
library(dplyr)
library(SeqArray)
library(SNPRelate)

hap.fn <- snakemake@input[["hapf"]]
leg.fn <- snakemake@input[["legf"]]
map.f <- snakemake@input[["mapf"]]
sample.fn <- snakemake@input[["filt_f"]]
out.gdsfn <- snakemake@output[["gdsf"]]
out.geno_gdsfn <- snakemake@output[["geno_gdsf"]]
chrom <- snakemake@params[["chrom"]]
compress.geno="LZ4_RA.fast"
compress.annotation="LZ4_RA.fast"






map_df <- read_delim(map.f,col_names = c("SNP", "pos", "map"),
                     delim = " ") %>%
    mutate(chrom = as.character(chrom))
sample.id <- scan(sample.fn,what=character(),sep="\n")
sample.id_haplo <- c(sapply(sample.id,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))
leg_df <- read_delim(leg.fn,delim=" ") %>%
    mutate(allele=paste0(allele0,",",allele1),chrom=as.character(chrom)) %>%
    mutate(snp_id=1:n()) %>% rename(SNP=ID)
leg_df <- inner_join(map_df,leg_df)

haplo_mat <- read_delim(hap.fn,delim=" ",col_names = F)
haplo_mat <- data.matrix(haplo_mat)[leg_df$snp_id,]

geno_mat <- haplo_2_geno_i(haplo_mat,snps_in_rows = T)
p <- nrow(geno_mat)
stopifnot(ncol(geno_mat)==length(sample.id))
stopifnot(ncol(haplo_mat)==length(sample.id_haplo))
gc()
tfile <- tempfile()
tfile2 <- tempfile()
snpgdsCreateGeno(gds.fn = tfile,genmat = haplo_mat,
                 sample.id = sample.id_haplo,
                 snp.id = 1:p,
                 snp.rs.id = leg_df$SNP,
                 snp.chromosome = leg_df$chrom,
                 snp.position = leg_df$pos,
                 snp.allele = leg_df$allele,
                 compress.annotation = compress.annotation,
                 compress.geno = compress.geno)

snpgdsCreateGeno(gds.fn = tfile2,genmat = geno_mat,
                 sample.id = sample.id,
                 snp.id = 1:p,
                 snp.rs.id = leg_df$SNP,
                 snp.chromosome = leg_df$chrom,
                 snp.position = leg_df$pos,
                 snp.allele = leg_df$allele,
                 compress.annotation = compress.annotation,
                 compress.geno = compress.geno)


seqSNP2GDS(gds.fn =tfile ,out.fn =out.gdsfn ,storage.option ="LZ4_RA.fast")
seqSNP2GDS(gds.fn =tfile2 ,out.fn =out.geno_gdsfn ,storage.option ="LZ4_RA.fast")
# gds <- seqOpen(out.gdsfn,readonly=F)
# gdsfmt::add.gdsn(index.gdsn(gds,"annotation/info"),"map",leg_df$map,replace=T,compress="LZ4_RA.fast")
# seqClose(gds)
