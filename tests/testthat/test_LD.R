context("Test LD calculations")

test_that("theta is computed correctly",{
  m <- 100
  nmsum = sum(1 / (1:(2*m-1)))
  theta = (1/nmsum) / (2*m + 1/nmsum)
  expect_equal(calc_theta(m),theta)
})





test_that("Equal to R implementation",{
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  mapdat <- cumsum(runif(p))

  tLD <- LDshrink(haplomat,mapdat,na.rm = F)
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


test_that("na.rm=T == na.rm==F",{
  # library(LDshrink)
  n <- 500
  p <- 1100

  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  mapdat <- cumsum(runif(p))
  #data("haplomat")
  
  #data("mapdat")
  ahaplomat <- haplomat+0
  bhaplomat <- haplomat+0
  tLD <- LDshrink(ahaplomat,mapdat)
  expect_equal(ahaplomat,haplomat)
  fLD <- LDshrink(bhaplomat,mapdat,na.rm=F)
  expect_equal(bhaplomat,haplomat)
  aLD <- LDshrink(ahaplomat,mapdat,m,Ne,cutoff,na.rm=F,useAlt = T)
  tmb <- microbenchmark::microbenchmark(aLD=LDshrink(bhaplomat,mapdat,na.rm=F,useAlt=T),times=50)
  expect_equal(tLD,fLD)
  expect_equal(aLD,tLD)
})


test_that("alternative shrinkage",{
  p <- 4
  mapdat <- cumsum(runif(p))
  m=85
  Ne=11490.672741
  cutoff = 0.001
  adist <- altDist(mapdat,m,Ne,cutoff)
  tdist <- trueDist(mapdat,m,Ne,cutoff)
  expect_equal(adist,tdist)
  
  
  
})

test_that("Equal to R implementation",{
  n <- 1000
  p <- 50
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  mapdat <- cumsum(runif(p))
  #data("haplomat")
  
  #data("mapdat")
  
  tLD <- LDshrink(haplomat,mapdat,na.rm=T)
  

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











