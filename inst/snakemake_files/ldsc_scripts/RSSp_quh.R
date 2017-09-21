library(RSSp)
library(readr)
library(dplyr)


quhf <- snakemake@input[["quhf"]]
traitn <- snakemake@params[["traitn"]]
stopifnot(file.exists(quhf))
quh_df <- read_delim(quhf,delim="\t")

p <- nrow(quh_df)
n <- unique(quh_df$N)
p_n <- p/n
rssp_res <- RSSp(fgeneid=traitn, D=quh_df$D, quh=quh_df$quh,p_n=p_n)
write_delim(rssp_res,snakemake@output[["RSSp_resf"]],delim="\t")


