context("Test LD calculations")
test_that("theta is computed correctly", {
  m <- 100
  nmsum <- sum(1 / (1:(2*m-1)))
  theta <- (1/nmsum) / (2*m + 1/nmsum)
  expect_equal(calc_theta(m), theta)
})

test_that("Dense ldshrink is equal to a toy R implementation", {
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
  haplo_bmat <- matrix(bit::as.bit(haplomat),nrow = n,ncol = p)
  mapdat <- cumsum(runif(p))

  tLD <- ldshrink(haplomat, mapdat, na.rm = F)
  calcLDR <- function(hmata, mapa, m=85, Ne=11490.672741, cutoff = 0.001,isgeno = FALSE){
    S <- stats::cov(hmata)
    p <- length(mapa)
    td <- abs(outer(mapa, mapa, `-`))
     if(isgeno){
      S <- 0.5*S
    }
  
    rho  <-  4*Ne*(td)/100;
    rho <- -rho/(2*m);
    tshrinkage <- exp(rho);
    tshrinkage[tshrinkage<cutoff] <- 0
    # diag(tshrinkage) <- 1
    S <- S*tshrinkage
    theta <- calc_theta(m)
    
    eye <- diag(p)*(0.5*theta * (1-0.5*theta))
    SigHat <-  ((1-theta)*(1-theta))*S+eye
    return(stats::cov2cor(SigHat))
  }
  
  RLD <- calcLDR(haplomat, mapdat)
  expect_equal(RLD, tLD)
})




test_that("Equal to R implementation", {
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)
  
  tLD <- ldshrink(haplomat,mapdat,na.rm = F)
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


test_that("Equal to R implementation", {
  n <- 1000
  p <- 50
  
  haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
  mapdat <- cumsum(runif(p))
  #data("haplomat")
  
  #data("mapdat")
  
  tLD <- ldshrink(haplomat,mapdat,na.rm=T)
  

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











