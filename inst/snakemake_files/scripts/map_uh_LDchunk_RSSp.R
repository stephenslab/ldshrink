library(RcppEigenH5)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
LDchunk <- as.integer(as.numeric(snakemake@params[["LDchunk"]]))
LDchunkf <- snakemake@output[["LDchunkf"]]
ymatf <- snakemake@input[["ymatf"]]
cores <- snakemake@threads




gds <- seqOpen(gdsf, readonly = T)

LD_region <- seqGetData(gds, "annotation/info/LD_chunk")
good_LD <- LD_region == LDchunk
stopifnot(sum(good_LD) > 0)
seqSetFilter(gds, variant.sel = good_LD)
snp_id <- seqGetData(gds,var.name="variant.id")
tparam_df <- read_df_h5(ymatf, "SimulationInfo")
ymat <- read_2d_mat_h5(ymatf, "trait", "ymat")
#p <- length(LD_region)

uh_mat <- gen_bhat_se_block_gds(gds,
                                ymat,
                                cores,
                                tparam_df, na.rm = T)

write_mat_h5(LDchunkf, as.character(LDchunk),
             "uh",
             data = uh_mat,
             deflate_level = 0,
             doTranspose = F)
#snp_id <- as.integer((1:p))[good_LD]
write_ivec_h5(LDchunkf,
              as.character(LDchunk),"snp_id",data = snp_id,deflate_level = 4L)
write_df_h5(tparam_df,"SimulationInfo",LDchunkf)
