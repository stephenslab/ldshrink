# Download `genotype2.mat` from the following link:
#https://uchicago.app.box.com/v/example2

library(RcppEigenH5)
shrink_R <- read_2d_mat_h5("/home/nwknoblauch/Downloads/genotype2.mat","/","shrink_R")
devtools::use_data(shrink_R)


panel_mapfile <- "/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/interpolated_from_hapmap/chr19.interpolated_genetic_map.gz"
map_dat <- readr::read_delim(panel_mapfile,delim=" ",col_names=c("rsid","pos","map"))
devtools::use_data(map_dat)

liftover_allf <- "~/Dropbox/ldetect-data/EUR/fourier_ls-all.bed"
break_df <- readr::read_delim(liftover_allf,delim="\t",trim_ws = T)

liftover_bed <- mutate(liftover_bed,bp=as.numeric(stop-start),ldp=(bp/1000)^2)


all_ldmats <- sapply(liftover_bed$bp,function(x){matrix(0,x/1000,x/1000)},simplify = F)

liftover_coord <- mutate(liftover_bed,chr=as.integer(gsub("chr","",chr))) %>% gather(start_stop,value,-chr) %>% arrange(chr,value) %>% mutate(pos=1:n())
ggplot(liftover_coord,aes(x=pos,y=value,col=factor(chr)))+geom_point()

test_vcf <- "~/Dropbox/LDshrink/inst/test_gdsf/sub_19.recode.vcf"
test_gds <- "~/Dropbox/LDshrink/inst/test_gdsf/sub_19.gds"

test_hapf <- "~/Dropbox/LDshrink/inst/test_gdsf/sub_19.impute.hap"
test_legf <- "~/Dropbox/LDshrink/inst/test_gdsf/sub_19.impute.legend"
test_sampf <- "~/Dropbox/LDshrink/inst/test_gdsf/sub_19.impute.hap.indv"

sample.id <- scan(test_sampf,what=character(),sep="\n")
sample.id_haplo <- c(sapply(sample.id,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))
chrom <- "19"
leg_df <- read_delim(test_legf,delim=" ") %>%
  mutate(allele=paste0(allele0,",",allele1),chrom=as.character(chrom)) %>%
  mutate(snp_id=1:n()) %>% rename(SNP=ID)

haplo_mat <- data.matrix(read_delim(test_hapf,delim=" ",col_names = F,col_types = cols(.default = col_double())))

tt_geno_mat <- apply(haplo_mat,1,function(x){2-sapply(split(x,gl(n = length(x)/2,k = 2,length = length(x))),sum)})
geno_mat <- t(2-haplo_2_geno(haplo_mat,snps_in_rows = T))
max()
ogds <- seqVCF2GDS(test_vcf,test_gds,optimize = F)
gds <- seqOpen(test_gds)
t_dosage_mat <- seqGetData(gds,"$dosage")
max(t_dosage_mat-geno_mat)
p <- nrow(geno_mat)
stopifnot(ncol(geno_mat)==length(sample.id))
gc()
tfile <- tempfile()
snpgdsCreateGeno(gds.fn = tfile,genmat = geno_mat,
                 sample.id = sample.id,
                 snp.id = 1:p,
                 snp.rs.id = leg_df$SNP,
                 snp.chromosome = leg_df$chrom,
                 snp.position = leg_df$pos,
                 snp.allele = leg_df$allele,
                 compress.annotation = compress.annotation,
                 compress.geno = compress.geno)

seqSNP2GDS(gds.fn =tfile ,out.fn =out.gdsfn ,storage.option ="LZ4_RA.fast")
