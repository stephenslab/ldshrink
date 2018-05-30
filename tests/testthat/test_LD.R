context("Test LD calculations")

test_that("theta is computed correctly",{
  m <- 100
  nmsum = sum(1 / (1:(2*m-1)))
  theta = (1/nmsum) / (2*m + 1/nmsum)
  expect_equal(calc_theta(m),theta)
})





test_that("Equal to R implementation",{
  data("haplomat")
  data("mapdat")

  tLD <- LDshrink(haplomat,mapdat)
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
  
  RLD <- calcLDR(haplomat,mapdat)
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









