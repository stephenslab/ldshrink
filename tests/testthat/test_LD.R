context("Test LD calculations")



# 
# 
# test_that("HDF5 backed data",{
#   
#   
#   library(LDshrink)
#   library(tidyverse)
#   library(testthat)
#   library(SeqSupport)
#   m=85
#   Ne=1490.672741
#   cutoff=1e-3
#   data("haplomat")
#   data("mapdat")
#   Hpanel <- haplomat
#   tmap <- mapdat
# 
#   
#   p <- ncol(Hpanel)
#   
#   
#   tfa <- tempfile()
#   tfo <- tempfile()
#   EigenH5::write_matrix_h5(tfa,"/","dosage",t(Hpanel))
#   EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
#   n <- nrow(Hpanel)
#   tby <- as.integer(100)
#   snp_start <- as.integer(seq(0,p-1,by=tby))
#   snp_stop <- as.integer(pmin(snp_start+tby,p))
#   snp_size <- as.integer(snp_stop-snp_start)
#   
#   snp_dff <- data_frame(filenames=tfa,
#                         groupnames="/",
#                         datanames="dosage",
#                         row_offsets=snp_start,
#                         row_chunksizes=snp_size,
#                         col_offsets=0L,
#                         col_chunksizes=as.integer(n)) %>% mutate(chunk_group=1:n())
#   map_dff <- mutate(snp_dff,col_chunksizes=1L,groupnames="SNPinfo",datanames="map")
#   data_dff <- bind_rows(snp_dff,map_dff)
# 
#   sdfl <- list()
#   sdfl[["Q"]] <- match_df_h5(snp_dff,
#                        output_filenames=tfo,
#                        output_group_prefix="EVD",
#                        output_datanames="Q")
#   sdfl[["D"]] <- match_df_h5(snp_dff,
#                        output_filenames=tfo,
#                        output_group_prefix="EVD",
#                        output_datanames="D",isVector=T)
#   sdfl[["R"]] <- match_df_h5(snp_dff,
#                              output_filenames=tfo,
#                              output_group_prefix="LD",
#                              output_datanames="R")
#   sdfl[["L2"]] <- match_df_h5(snp_dff,
#                              output_filenames=tfo,
#                              output_group_prefix="L2",
#                              output_datanames="L2",isVector=T)
# 
#   out_df <- bind_rows(sdfl)
#   dplyr::filter(out_df,!create_dynamic) %>% EigenH5::create_mat_l()
#   
#   calc_LD_chunk_h5(input_dff=data_dff,
#                    output_dff=out_df,
#                    m=m,Ne=Ne,
#                    cutoff=cutoff,
#                    SNPfirst=T,
#                    evd=T,
#                    svd=F,
#                    df=F,
#                    r2cutoff=0.01)
#   
#   
#   
#  
#   nS <- EigenH5::read_mat_h5(tfo,"LD/1","R")
#   nRsig_1 <- calcLD(Hpanel[,1:100],tmap[1:100],m,Ne,cutoff)
#   expect_equal(nS,nRsig_1)
#   evdR <- eigen(nRsig_1)
#   nD <- EigenH5::read_vector_h5(tfo,"EVD/1","D")
#   expect_equal(evdR$values,nD)
#   
#   nQ <- EigenH5::read_mat_h5(tfo,"EVD/1","Q")
#   nnS <- nQ%*%diag(nD)%*%t(nQ)
#   expect_equal(nnS,nS)  
#   Revd <- evdR$vectors%*%diag(evdR$values)%*%t(evdR$vectors)
#   expect_equal(nRsig_1,Revd)
#   for(i in 1:100){
#     expect_equal(c(nQ[,i]%*%nQ[,i]),1)
#     if(i>1){
#       expect_equal(c(nQ[,i-1]%*%nQ[,i]),0)
#     }
#   }
# })
# 
# 
# test_that("HDF5 backed data with one chunk",{
#   
#   
#   library(LDshrink)
#   library(EigenH5)
#   m=85
#   Ne=1490.672741
#   cutoff=1e-3
#   data("haplomat")
#   data("mapdat")
#   Hpanel <- haplomat
#   tmap <- mapdat
#   
#   
#   p <- ncol(Hpanel)
#   
#   
#   tfa <- tempfile()
#   tfo <- tempfile()
#   EigenH5::write_matrix_h5(tfa,"/","dosage",Hpanel)
#   #EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
#   # ld_reg <- c(rep(1,p/2),rep(2,p/2))
#   EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
#   dosage_input_dff <- tibble::data_frame(filenames=tfa,groupnames="/",datanames="dosage",chunk_group=1L)
#   map_input_dff <- tibble::data_frame(filenames=tfa,groupnames="SNPinfo",datanames="map",chunk_group=1L)
#   input_ff <- rbind(map_input_dff,dosage_input_dff)
#   
#   R_ff <- tibble::data_frame(
#     filenames=c(tfo),
#     groupnames=c("R/1","Q/1","D/1","L2/1"),
#     datanames=c("R","Q","D","L2"),
#     row_chunksizes=c(p,p,p,p),
#     col_chunksizes=c(p,p,1,1),
#     row_c_chunksizes=c(p,p,1,1),
#     col_c_chunksizes=c(p,p,1,1),datatypes="numeric",chunk_group=1L)
#   create_mat_l(R_ff)
#   
#   
#   calc_LD_chunk_h5(input_ff,R_ff,m=m,Ne=Ne,cutoff=cutoff,SNPfirst=F)
#   nS <- EigenH5::read_mat_h5(tfo,"R/1","R")
#   rR <- calcLD(Hpanel,tmap,m,Ne,cutoff)
#   testthat::expect_equal(nS,rR)
#  
#   evdR <- eigen(nS)
#   nD <- EigenH5::read_vector_h5(tfo,"D/1","D")
#   expect_equal(evdR$values,nD)
#   
#   nQ <- EigenH5::read_mat_h5(tfo,"Q/1","Q")
#   nnS <- nQ%*%diag(nD)%*%t(nQ)
#   expect_equal(nnS,nS)  
# 
# })






test_that("Equal to R implementation",{
  
  

  tLD <- LDshrink(Hpanel,mapdat)
  calcLDR <- function(hmata,mapa,m=85,Ne=11490.672741, cutoff = 0.001){
    S <- cov(hmata)
    p <- length(mapa)
    td <- abs(outer(mapa,mapa,`-`))
    # td[lower.tri(td)] <- 0
    # td <- td+t(td)
    rho = 4*Ne*(td)/100;
    rho=-rho/(2*m);
    tshrinkage=exp(rho);
    tshrinkage[tshrinkage<cutoff] <- 0
    diag(tshrinkage) <- 1
    S <- S*tshrinkage
    theta <- calc_theta(m)
    
    eye <- diag(p)*(0.5*theta * (1-0.5*theta))
    SigHat <-  ((1-theta)*(1-theta))*S+eye
    return(cov2cor(SigHat))
  }
  
  RLD <- calcLDR(Hpanel,mapdat)
  expect_equal(RLD,tLD)  
})

test_that("blockwise implementation works",{
  mdir <- system.file("test_gdsf/test_SNP.h5",package="LDshrink")
  of <- tempfile()
  snp_df <- EigenH5::read_df_h5(h5filepath = mdir,groupname = "SNPinfo")
  p <- nrow(snp_df)
  snp_df <- mutate(snp_df,region_id=sort(rep(1:2,p/2)))
  
  tr <- chunkwise_LD_h5(mdir,of,snp_df)
  ls_h5(of,"LD/1")
  rR <- read_matrix_h5(of,"LD/1","R")
  
})









