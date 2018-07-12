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




test_that("Sparse and Dense implementations are equivalent for non-LDshrink",{
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)
  
  tLD <- cor(haplomat)
  sLD <- sparse_LDshrink(scaled_data = scale(haplomat,center=T,scale=F),
                         mapd = mapdat,m=85,Ne=11490.672741, cutoff = 0.001,useLDshrink = F)
  dsLD <- as.matrix(sLD)
  
  
  # mb <- microbenchmark::microbenchmark(tLD = LDshrink(haplomat,mapdat,na.rm = F),
  #                                      sLD = sparse_LDshrink(scaled_data = scale(haplomat,center=T,scale=F),
  #                                                             mapd = mapdat,m=85,Ne=11490.672741, cutoff = 0.001))
  expect_equivalent(tLD,dsLD)
})




test_that("Sparse and Dense implementations are equivalent",{
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)
  
  tLD <- LDshrink(haplomat,mapdat,na.rm = F)
  sLD <- sparse_LDshrink(scaled_data = scale(haplomat,center=T,scale=F),
                         mapd = mapdat,m=85,Ne=11490.672741, cutoff = 0.001)
  dsLD <- as.matrix(sLD)
  
  
  # mb <- microbenchmark::microbenchmark(tLD = LDshrink(haplomat,mapdat,na.rm = F),
  #                                      sLD = sparse_LDshrink(scaled_data = scale(haplomat,center=T,scale=F),
  #                                                             mapd = mapdat,m=85,Ne=11490.672741, cutoff = 0.001))
  expect_equivalent(tLD,dsLD,1e-3)
})



test_that("Equal to R implementation",{
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)
  
  tLD <- LDshrink(haplomat,mapdat,na.rm = F)
  calcLDR <- function(hmata,mapa,m=85,Ne=11490.672741, cutoff = 0.001){
    S <- stats::cov(hmata)
    p <- length(mapa)
    # nmap <- -((4*Ne*mapa)/100)/(2*m)
    # ntd <- outer(nmap,nmap,function(x,y)exp(x-y))
    # nntd <- outer(snmap,snmap,function(x,y)ifelse(x>y,y/x,x/y))
    # ntd[is.na(ntd)] <- 0
    # diag(ntd) <- 1
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


test_that("LD2df correctly subsets", {
  n <- 500
  p <- 1100
  p_tot <- (p^2-p)/2+p
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)

  tLDdf <- ld2df(LDshrink(haplomat,mapdat,na.rm = F,m = 85,Ne=11490.672741,cutoff=0.001),rsid = as.character(1:p),r2cutoff = 0.001)
  nhaplomat <- scale(haplomat,center=T,scale=F)
  nLDdf <- ld2df_p(scaled_data = nhaplomat,mapd = mapdat,rsid =as.character(1:p),m=85,Ne=11490.672741, cutoff = 0.001,r2cutoff = 0.001)
  
  nR <- matrix(NA,p,p)
  nR[cbind(as.integer(tLDdf$rowsnp),as.integer(tLDdf$colsnp))] <- tLDdf$r2
  cR <- matrix(NA,p,p)
  cR[cbind(as.integer(nLDdf$rowsnp),as.integer(nLDdf$colsnp))] <- nLDdf$r2
  testthat::expect_equal(nrow(nLDdf),nrow(tLDdf))
  testthat::expect_equal(cR[upper.tri(cR)],nR[upper.tri(nR)],1e-4)
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











