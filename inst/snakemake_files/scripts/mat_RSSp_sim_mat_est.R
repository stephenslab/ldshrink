library(RSSp)
library(dplyr)
library(purrr)
library(readr)
rds_data <- readRDS(snakemake@input[["rdsf"]])
outf <- snakemake@output[["dff"]]

quh_mat <- rds_data[["quh_mat"]]
D <- rds_data[["D"]]
tparam_df <- rds_data[["tparam_df"]]

colnames(quh_mat) <- as.character(tparam_df$fgeneid)


n <-unique(tparam_df$n)
stopifnot(length(n)==1)

rss_res <- cross(list(doConfound=c(T,F),log_params=c(T,F),useGradient=c(T))) %>%  invoke_map_dfr("RSSp_run_mat_quh",.,quh_mat_d=quh_mat,D=D,n=n) %>% inner_join(tparam_df)

write_delim(rss_res,outf,delim="\t")
