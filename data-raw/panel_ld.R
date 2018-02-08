# Download `genotype2.mat` from the following link:
#https://uchicago.app.box.com/v/example2

library(RcppEigenH5)
library(tidyverse)
shrink_R <- read_2d_mat_h5("/home/nwknoblauch/Downloads/genotype2.mat","/","shrink_R")
devtools::use_data(shrink_R)


panel_mapfile <- "/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/interpolated_from_hapmap/chr19.interpolated_genetic_map.gz"
map_dat <- readr::read_delim(panel_mapfile,delim=" ",col_names=c("rsid","pos","map"))
devtools::use_data(map_dat)

liftover_allf <- "~/Dropbox/ldetect-data/EUR/fourier_ls-all.bed"
break_df <- read_delim(liftover_allf,delim="\t",trim_ws = T)
break_df <- dplyr::mutate(break_df,chr=as.integer(gsub("chr","",chr))) %>% mutate(region_id=1:n())
devtools::use_data(break_df)
liftover_bed <- mutate(liftover_bed,bp=as.numeric(stop-start),ldp=(bp/1000)^2)


# all_ldmats <- sapply(liftover_bed$bp,function(x){matrix(0,x/1000,x/1000)},simplify = F)

liftover_coord <- mutate(liftover_bed,chr=as.integer(gsub("chr","",chr))) %>% gather(start_stop,value,-chr) %>% arrange(chr,value) %>% mutate(pos=1:n())
ggplot(liftover_coord,aes(x=pos,y=value,col=factor(chr)))+geom_point()

