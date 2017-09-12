
library(dplyr)
library(RcppEigenH5)
library(purrr)
library(progress)

inf <- snakemake@input[["evdf"]]
outf <- snakemake@output[["evd_chunkf"]]

nchunks <- length(inf)
pb <- progress_bar$new(total = nchunks)
ld_l <- list()
mat_l <- list()
D_l <- list()
save.image()
for (i in 1:nchunks){
    x <- inf[i]
    ldi <- read_df_h5(x, "LDinfo")
    D <-  read_dvec(x, "EVD", "D")
    Q <- read_2d_mat_h5(x, "EVD", "Q")

    region_id <- as.character(unique(ldi$region_id))
    write_mat_h5(outf, region_id, "Q", Q, 0)
    write_dvec_h5(outf, region_id, "D", D, 0)
    ld_l[[i]] <- ldi
    pb$tick()
}
    
ld_info <- bind_rows(ld_l)

write_df_h5(ld_info, "LDinfo", outf)
