#download omni genetic map
ftp_omni <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/"
ftp_hapmap <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/"

populations <- c("ACB","ASW","CDX","CEU","CHB","CHS","CLM","FIN","GBR","GIH","IBS","JPT","KHV","LWK","MKK","MXL","PEL","PUR","REA","TSI","YRI")



pop_size_url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/README_omni_recombination_20130507"
ind_file <- "ftp:/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/omni_haplotypes/20111117_omni_sample_list.txt"
n_ind_file <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130501_sample_info/20130501_sample_info.txt"
ind_d <- scan(ind_file,what=character(),sep="\n")
pop_size_df <- readr::read_delim(pop_size_url,skip = 24,col_names=c("Pop","Ne"),delim="\t") 
pop_df <- read_delim(n_ind_file,delim="\t")
colnames(pop_df) <- gsub(pattern = " ","_",fixed=T,colnames(pop_df))
p1_df <- filter(pop_df,Phase1==1)

p1_popsize <- filter(pop_df,Phase1==1) %>% group_by(Population_Code) %>% summarise(m=n()) %>% rename(Pop=Population_Code)
map_parameters <- inner_join(p1_popsize,pop_size_df)

devtools::use_data(map_parameters,compress = "gzip")  

# LDshrink <- function(haplo_panel,map_population="CEU",map_data=map_fun(map_population),m=m_fun(map_population),Ne=Ne_fun(map_population),cutoff=1e-3,isGeno=NA,cov_2_cor=T,na.rm=T){
# X <- read.table("some_file.txt")
# R <- LDshrink(X)
# R <- LDshrink(X,"YRI")
# R <- LDshrink(X,"Other_Pop",map_data=map,m=m,Ne=Ne)
library(curl)
library(tidyverse)
library(stringr)
# url = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/11/PXD000299/"
h <-  new_handle(dirlistonly = TRUE)
con <-  curl(ftp_hapmap, "r", h)
tbl_hapmap <-  read.table(con, stringsAsFactors = TRUE, fill = TRUE,col.names = "filename")
hapmap_df <- tbl_hapmap %>%
  filter(str_detect(filename,"gz")) %>%
  mutate(file_url=paste0(ftp_hapmap,filename))
curl_fetch_disk(url = hapmap_df$file_url,path = paste0("inst/genetic_map/",hapmap_df$filename))
close(con)

con <-  curl(ftp_omni, "r", h)
tbl_omni <- read.table(con, stringsAsFactors = TRUE, fill = TRUE)
close(con)


urls_ <- paste0(url, tbl[1:5,1])
fls = basename(urls)
curl_fetch_disk(urls[1], fls[1])


parse_map <- function(panel_name="CEU"){
  inst_dir <- "~/Dropbox/projectile/LDshrink/inst/genetic_map/CEU"
  
  map_file_df <- data_frame(filename=dir(inst_dir,full.names = T)) %>% 
    mutate(chrom=as.integer(gsub(paste0(panel_name,"-([0-9]+)-final.txt.gz"),"\\1",basename(filename))))
  CEU_map <- pmap_dfr(map_file_df,function(filename,chrom){
  read_delim(filename,delim="\t",col_names=c("pos","rate","map","filtered"),skip = 1,trim_ws = T) %>% mutate(chrom=chrom)
  }) %>% arrange(chrom,pos)
  
devtools::use_data(CEU_map,compress = "gzip")  
  
}

