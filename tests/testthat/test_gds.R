context("gds files")

test_that("Getting from haplotype to genotype works as expected",{
  inst_dir <- system.file("test_gdsf",package="LDshrink")
  test_vcf <- file.path(inst_dir,"sub_19.recode.vcf")
  test_gds <- file.path(inst_dir,"sub_19.gds")
  test_mapf <-file.path(inst_dir,"map_sub_19.txt.gz") 
  test_hapf <- file.path(inst_dir,"sub_19.impute.hap")
  test_legf <- file.path(inst_dir,"sub_19.impute.legend")
  test_sampf <- file.path(inst_dir,"sub_19.impute.hap.indv")
  
  sample.id <- scan(test_sampf,what=character(),sep="\n")
  sample.id_haplo <- c(sapply(sample.id,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))
  chrom <- "19"
  leg_df <- readr::read_delim(test_legf,delim=" ") %>%
    mutate(allele=paste0(allele0,",",allele1),chrom=as.character(chrom)) %>%
    mutate(snp_id=1:n()) %>% rename(SNP=ID)
  map_df <- read_delim(test_mapf,delim=" ",col_names = c("SNP","pos","map")) %>% mutate(chrom="19")
  
  haplo_mat <- data.matrix(read_delim(test_hapf,delim=" ",col_names = F,col_types = cols(.default = col_double())))
  
  tt_geno_mat <- apply(haplo_mat,1,function(x){2-sapply(split(x,gl(n = length(x)/2,k = 2,length = length(x))),sum)})
  dimnames(tt_geno_mat)<- NULL
  geno_mat <- t(2-haplo_2_geno(haplo_mat,snps_in_rows = T))
  ogeno_mat <- haplo_2_geno(haplo_mat,snps_in_rows = T)
  expect_equal(geno_mat,tt_geno_mat)
  
  
  ogds <- seqVCF2GDS(test_vcf,test_gds,optimize = F)
  gds <- seqOpen(test_gds)
  t_dosage_mat <- seqGetData(gds,"$dosage")
  dimnames(t_dosage_mat) <- NULL
  expect_equal(geno_mat,t_dosage_mat)

  p <- nrow(ogeno_mat)
  stopifnot(ncol(ogeno_mat)==length(sample.id))
  gc()
  tfile <- tempfile()
  snpgdsCreateGeno(gds.fn = tfile,genmat = ogeno_mat,
                   sample.id = sample.id,
                   snp.id = 1:p,
                   snp.rs.id = leg_df$SNP,
                   snp.chromosome = leg_df$chrom,
                   snp.position = leg_df$pos,
                   snp.allele = leg_df$allele)
  otfile <- tempfile()
  seqSNP2GDS(gds.fn =tfile ,out.fn =otfile ,storage.option ="LZ4_RA.fast")
  ngds <- seqOpen(otfile)
  ndosage <-   seqGetData(ngds,"$dosage")
  expect_equal(ndosage,seqGetData(gds,"$dosage"))
  
  
})