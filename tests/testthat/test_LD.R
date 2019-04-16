context("Test LD calculations")
# test_that("theta is computed correctly", {
# 
#   expect_equal(calc_theta(m), theta)
# })

test_that("ldshrink can work like base R for sample correlation",{
  # n <- 500
                                        # p <- 1100
    library(ldshrink)
  data("reference_genotype")
  data("reference_map")
  
  rResult <- estimate_LD(reference_panel = reference_genotype[,1:5],method="sample",map = reference_map[1:5],output = "matrix")
  
  R <- cor(reference_genotype[,1:5])
  expect_equal(R,rResult,check.attrubites=F)

})

test_that("Dense ldshrink is equal to a toy R implementation", {

  n <- 500
  p <- 1100
  haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
  mapdat <- cumsum(runif(p))
  
  # haplomat <- matrix(sample(as.numeric(0:1), n*p, replace = T), n, p)
  # mapdat <- cumsum(runif(p))
  haplomat[,1:2] <- 1
  calcLDR <- function(hmata, mapa, m=85, Ne=11490.672741, cutoff = 0.001,isgeno=FALSE){
    S <- stats::cov(hmata)
    p <- length(mapa)
    td <- abs(outer(mapa, mapa, `-`))
    if(isgeno){
      S <- 0.5*S
    }
    # td[lower.tri(td)] <- 0
    # td <- td+t(td)
    rho  <-  4*Ne*(td)/100;
    rho <- -rho/(2*m);
    tshrinkage <- exp(rho);
    tshrinkage[tshrinkage<cutoff] <- 0
    #diag(tshrinkage) <- 1
    S <- S*tshrinkage
    theta <- ldshrink:::calc_theta(m)
    
    eye <- diag(p)*(0.5*theta * (1-0.5*theta))
    SigHat <-  ((1-theta)*(1-theta))*S+eye
   
    return(stats::cov2cor(SigHat))
  }
  RLD <- calcLDR(haplomat, mapdat)
  tLD <- estimate_LD(reference_panel = haplomat,map =  mapdat,output="matrix",isGenotype = F)
  expect_equal(RLD, tLD,check.attributes = F)
})




test_that("Dense ldshrink is equal to a toy R implementation for a relatively large  sample data", {

  N <- 2000
  p <- 2100
  x <- matrix(sample(0.0:2.0,size=N*p,replace=T),nrow = N)
  
  data("reference_map")
  n_map_pos <- sort(c(1,sample(2:(p-1),length(reference_map)-2,replace=F),p))
  n_map <- interpolate_genetic_map(reference_map,n_map_pos,1:p)
  expect_false(is.unsorted(n_map,strictly = T))
  ld_dist <- function(map_dist,m=85, Ne=11490.672741, cutoff = 0.001){
    rho  <-  4*Ne*(map_dist)/100;
    rho <- -rho/(2*m)
    tshrinkage <- exp(rho)
    tshrinkage[tshrinkage<cutoff] <- 0
    return(tshrinkage)
  }
  calcLDR <- function(hmata, mapa, m=85, Ne=11490.672741, cutoff = 0.001,isgeno=FALSE){
    S <- stats::cov(hmata)
    p <- length(mapa)
    td <- abs(outer(mapa, mapa, `-`))
    if(isgeno){
      S <- 0.5*S
    }
    tshrinkage <- ld_dist(td,m,Ne,cutoff)
    S <- S*tshrinkage
    theta <- ldshrink:::calc_theta(m)
    
    eye <- diag(p)*(0.5*theta * (1-0.5*theta))
    SigHat <-  ((1-theta)*(1-theta))*S+eye
    
    return(stats::cov2cor(SigHat))
  }
  RLD <- calcLDR(x, n_map,isgeno=TRUE)
  tLD <- estimate_LD(reference_panel = x,map =  n_map,output="matrix",isGenotype = T)
  expect_equal(RLD, tLD,check.attributes = F)
})














