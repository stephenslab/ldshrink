context("Test LD calculations")


test_that("Covariance works as expected",{
  
  inp_a <- matrix(runif(90*10),90,10)
  #inp_b <- matrix(runif(90*9),90,9)  
  
  R_cov <- cov(inp_a)
  
  ld_cov <- calc_cov(inp_a)
  expect_equal(ld_cov,R_cov)
})

test_that("cov2cor works as expected",{
  inp_a <- matrix(runif(90*10),90,10)
  #inp_b <- matrix(runif(90*9),90,9)  
  
  R_cov <- cov(inp_a)
  
  ld_cov <- calc_cov(inp_a)
  expect_equal(ld_cov,R_cov)
  R_cor <- cor(inp_a)
  m_cor <- cov_2_cor(ld_cov)
  expect_equal(m_cor,R_cor)  
  
})

# test_that("covariance from MKL works",{
#   
#   inp_a <- matrix(rnorm(9000*9,mean = 3,sd = 2),9000,9)
#   
#   R_cov <- cov(inp_a)
#   
#   
#   mkl_cov <- cov_mkl(inp_a)
#   S <- mkl_cov$S
#   expect_equal(mkl_cov$S,R_cov)
#   
# })

test_that("HDF5 backed data",{
  
  
  library(LDshrink)
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat

  
  p <- ncol(Hpanel)
  
  
  tfa <- tempfile()
  tfo <- tempfile()
  EigenH5::write_matrix_h5(tfa,"/","dosage",t(Hpanel))
  #EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
  ld_reg <- c(rep(1,p/2),rep(2,p/2))
  #EigenH5::write_vector_h5(tfa,"SNPinfo","region_id",ld_reg)
  calc_ld_h5_exp(tfa,tfo,ld_reg,tmap,m = m,Ne=Ne,cutoff=cutoff,SNPfirst = T)
  nS <- EigenH5::read_mat_h5(tfo,"LD/1","R")
  nRsig_1 <- calcLD_prel(Hpanel[,1:(p/2)],tmap[1:(p/2)],m,Ne,cutoff)
  testthat::expect_equal(nS,nRsig_1)
  
  
  tfa <- tempfile()
  tfo <- tempfile()
  EigenH5::write_matrix_h5(tfa,"/","dosage",Hpanel)

  EigenH5::write_vector_h5(tfa,"SNPinfo","region_id",ld_reg)
  calc_ld_h5_exp(tfa,tfo,ld_region = ld_reg,mapd = tmap,m = m,Ne=Ne,cutoff=cutoff,SNPfirst = F)
  nS <- EigenH5::read_mat_h5(tfo,"LD/1","R")
  nRsig_1 <- calcLD_prel(Hpanel[,1:(p/2)],tmap[1:(p/2)],m,Ne,cutoff)
  expect_equal(nS,nRsig_1)
  evdR <- eigen(nRsig_1)
  nD <- EigenH5::read_vector_h5(tfo,"EVD/1","D")
  expect_equal(evdR$values,nD)
  
  nQ <- EigenH5::read_mat_h5(tfo,"EVD/1","Q")
  nnS <- nQ%*%diag(nD)%*%t(nQ)
  expect_equal(nnS,nS)  
  Revd <- evdR$vectors%*%diag(evdR$values)%*%t(evdR$vectors)
  expect_equal(nRsig_1,Revd)
  nS <- EigenH5::read_mat_h5(tfo,"LD/2","R")
  nRsig_2 <- calcLD_prel(Hpanel[,-(1:(p/2))],tmap[-(1:(p/2))],m,Ne,cutoff)
  testthat::expect_equal(nS,nRsig_2)
  
  evdR <- eigen(nRsig_2)
  nD <- EigenH5::read_vector_h5(tfo,"EVD/2","D")
  expect_equal(evdR$values,nD)
  
  nQ <- EigenH5::read_mat_h5(tfo,"EVD/2","Q")
  nnS <- nQ%*%diag(nD)%*%t(nQ)
  expect_equal(nnS,nS)  
  Revd <- evdR$vectors%*%diag(evdR$values)%*%t(evdR$vectors)
  expect_equal(nRsig_2,Revd)
})


test_that("HDF5 backed data with one chunk",{
  
  
  library(LDshrink)
  library(EigenH5)
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  p <- ncol(Hpanel)
  
  
  tfa <- tempfile()
  tfo <- tempfile()
  EigenH5::write_matrix_h5(tfa,"/","dosage",Hpanel)
  #EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
  # ld_reg <- c(rep(1,p/2),rep(2,p/2))
  EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
  dosage_input_dff <- tibble::data_frame(filenames=tfa,groupnames="/",datanames="dosage",chunk_group=1L)
  map_input_dff <- tibble::data_frame(filenames=tfa,groupnames="SNPinfo",datanames="map",chunk_group=1L)
  input_ff <- rbind(map_input_dff,dosage_input_dff)
  
  R_ff <- tibble::data_frame(
    filenames=c(tfo),
    groupnames=c("R/1","Q/1","D/1","L2/1"),
    datanames=c("R","Q","D","L2"),
    row_chunksizes=c(p,p,p,p),
    col_chunksizes=c(p,p,1,1),
    row_c_chunksizes=c(p,p,1,1),
    col_c_chunksizes=c(p,p,1,1),datatypes="numeric",chunk_group=1L)
  create_mat_l(R_ff)
  
  
  calc_LD_chunk_h5(input_ff,R_ff,m=m,Ne=Ne,cutoff=cutoff,SNPfirst=F)
  nS <- EigenH5::read_mat_h5(tfo,"R/1","R")
  rR <- calcLD(Hpanel,tmap,m,Ne,cutoff)
  testthat::expect_equal(nS,rR)
 
  evdR <- eigen(nS)
  nD <- EigenH5::read_vector_h5(tfo,"D/1","D")
  expect_equal(evdR$values,nD)
  
  nQ <- EigenH5::read_mat_h5(tfo,"Q/1","Q")
  nnS <- nQ%*%diag(nD)%*%t(nQ)
  expect_equal(nnS,nS)  

})



test_that("HDF5 backed data with 2 chunks",{
  
  
  library(LDshrink)
  library(EigenH5)
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  p <- ncol(Hpanel)
  
  
  tfa <- tempfile()
  tfo <- tempfile()
  EigenH5::write_matrix_h5(tfa,"/","dosage",Hpanel)
  #EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
  # ld_reg <- c(rep(1,p/2),rep(2,p/2))
  EigenH5::write_vector_h5(tfa,"SNPinfo","map",tmap)
  dosage_input_dff <- tibble::data_frame(filenames=tfa,groupnames="/",datanames="dosage",chunk_group=1L)
  map_input_dff <- tibble::data_frame(filenames=tfa,groupnames="SNPinfo",datanames="map",chunk_group=1L)
  input_ff <- rbind(map_input_dff,dosage_input_dff)
  
  R_ff <- tibble::data_frame(
    filenames=c(tfo),
    groupnames=c("R/1","Q/1","D/1","L2/1"),
    datanames=c("R","Q","D","L2"),
    row_chunksizes=c(p,p,p,p),
    col_chunksizes=c(p,p,1,1),
    row_c_chunksizes=c(p,p,1,1),
    col_c_chunksizes=c(p,p,1,1),datatypes="numeric",chunk_group=1L)
  create_mat_l(R_ff)
  
  
  calc_LD_chunk_h5(input_ff,R_ff,m=m,Ne=Ne,cutoff=cutoff,SNPfirst=F)
  nS <- EigenH5::read_mat_h5(tfo,"R/1","R")
  rR <- calcLD(Hpanel,tmap,m,Ne,cutoff)
  testthat::expect_equal(nS,rR)
  
  evdR <- eigen(nS)
  nD <- EigenH5::read_vector_h5(tfo,"D/1","D")
  expect_equal(evdR$values,nD)
  
  nQ <- EigenH5::read_mat_h5(tfo,"Q/1","Q")
  nnS <- nQ%*%diag(nD)%*%t(nQ)
  expect_equal(nnS,nS)  
  
})

# test_that("off-diagonal LD blocks can be calculated",{
#   
#   m=85
#   Ne=1490.672741
#   cutoff=1e-3
#   data("haplomat")
#   data("mapdat")
#   Hpanel <- haplomat
#   tmap <- mapdat
#   
#   
#   tRsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
#   theta <- calc_theta(m = m)
#   ldp <- c(m,Ne,cutoff,theta)
#   p <- ncol(Hpanel)
#   csize <- 501
#   chunkp <- as.integer(c(0L,1L,csize))
#   dRsig <- calcLD_par(hmat = Hpanel,map = tmap,ldparams = ldp,id = chunkp)
#   expect_true(all(abs(dRsig)<=1))
# 
#   # Rsig[lower.tri(Rsig)] <- 0
#   sub_tRsig <- tRsig[1:501,502:1000]
#   expect_equal(sub_tRsig,dRsig)
#   
# })


test_that("diagonal LD blocks can be calculated",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  tRsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  nRsig <- calcLD_prel(Hpanel,tmap,m,Ne,cutoff)
  expect_equal(nRsig,tRsig)
})

test_that("Cov2cor_p works as expected",{
  
  inp_a <- matrix(runif(90*10),90,10)
  inp_b <- matrix(runif(90*9),90,9)  
  
  R_cov <- cov(inp_a,inp_b)
  R_cor <- cor(inp_a,inp_b)
  
  
  ld_cov <- calc_cov_p(inp_a,inp_b)
  cova <- apply(inp_a,2,var)
  covb <- apply(inp_b,2,var)
  ld_cor <- cov_2_cor_exp_p(covmat = ld_cov,rowvar = cova,colvar = covb)
  expect_equal(ld_cor,R_cor)
  
})

test_that("Equal to R implementation",{
  
  
  
  tLD <- calcLDt(hmata = Hpanel,mapa = mapdat,cutoff = 0)
  tLD$SigHat[1:4,1:4]
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
  expect_equal(RLD,tLD$R)  
  
  
})









