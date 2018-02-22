# Generate test haplotype and genotype data (with and without missing data)




region_id <-  seqGetData(gds,"annotation/info/LD_chunk")
p <- nrow(test_gds_hap_leg)
test_gds_hap_ind <- scan(,what=character())
n <- length(test_gds_hap_ind)
hap_n <- n*2
 <- matrix(scan(file.path(test_gds_dir,"sub_19.impute.hap"),what = integer()),nrow = p,ncol = hap_n,byrow = T)

devtools::use_data()