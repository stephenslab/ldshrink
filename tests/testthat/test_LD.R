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
    S <- stats::cov(hmata)
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
    return(stats::cov2cor(SigHat))
  }
  
  RLD <- calcLDR(haplomat,mapdat)
  expect_equal(RLD,tLD)  
})



test_that("Equal to R implementation",{
  n <- 1000
  p <- 50
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  mapdat <- cumsum(runif(p))
  #data("haplomat")
  
  #data("mapdat")
  
  tLD <- LDshrink(haplomat,mapdat)
  
  LDshrink_alt <- function(haplo_panel,map_data,m=85,Ne=11490.672741,cutoff=1e-3,isGeno=NA,cov_2_cor=T,na.rm=T){
    if(is.na(isGeno)){
      isGeno <- max(haplo_panel,na.rm = na.rm)>1
    }
    stopifnot(!is.na(isGeno))
    Genomult <- ifelse(isGeno,0.5,1)
    haplo_panel <- scale(haplo_panel,center=T,scale=F)
    S <- stats::cov(haplo_panel,use=ifelse(na.rm,"complete.obs","all.obs"))*Genomult
    mS <- alt_shrinkCov(S,map_data,m,Ne,cutoff)
    if(cov_2_cor){
      return(stats::cov2cor(mS))
    }else{
      return(mS)
    }
  }
  nLD <- LDshrink_alt(haplomat,mapdat)
  
  expect_equal(tLD,nLD)
  
  calcLDR <- function(hmata,mapa,m=85,Ne=11490.672741, cutoff = 0.001){
    S <- stats::cov(hmata)
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
    return(stats::cov2cor(SigHat))
  }
  
  RLD <- calcLDR(haplomat,mapdat)
  expect_equal(RLD,tLD)  
})











