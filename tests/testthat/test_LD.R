context("Test LD calculations")
test_that("theta is computed correctly", {
  m <- 100
  nmsum <- sum(1 / (1:(2*m-1)))
  theta <- (1/nmsum) / (2*m + 1/nmsum)
  expect_equal(calc_theta(m), theta)
})

test_that("Equal to R implementation", {
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
  mapdat <- cumsum(runif(p))

  tLD <- LDshrink(haplomat, mapdat, na.rm = F)
  calcLDR <- function(hmata, mapa, m=85, Ne=11490.672741, cutoff = 0.001){
    S <- stats::cov(hmata)
    p <- length(mapa)
    td <- abs(outer(mapa, mapa, `-`))
    # td[lower.tri(td)] <- 0
    # td <- td+t(td)
    rho  <-  4*Ne*(td)/100;
    rho <- -rho/(2*m);
    tshrinkage <- exp(rho);
    tshrinkage[tshrinkage<cutoff] <- 0
    diag(tshrinkage) <- 1
    S <- S*tshrinkage
    theta <- calc_theta(m)
    
    eye <- diag(p)*(0.5*theta * (1-0.5*theta))
    SigHat <-  ((1-theta)*(1-theta))*S+eye
    return(stats::cov2cor(SigHat))
  }
  
  RLD <- calcLDR(haplomat, mapdat)
  expect_equal(RLD, tLD)
})




test_that("Sparse and Dense implementations are equivalent for non-LDshrink", {
  n <- 500
  p <- 1100
  
  haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)
  
  tLD <- cor(haplomat)
  sLD <- LDshrink::sparse_LDshrink(data = haplomat,
                                   mapd = mapdat,
                                   indices=0:(p-1),
                                   m=85,
                                   Ne=11490.672741,
                                   cutoff = 0.001,
                                   total_size=p,
                                   useLDshrink = F)
  dsLD <- as.matrix(sLD)
  

  expect_equivalent(tLD, dsLD)
})




test_that("Sparse and Dense implementations are equivalent", {
  n <- 5000
  p <- 1100
  
  haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
  nhaplomat <- haplomat+0
  tmap <- runif(p)
  
  mapdat <- cumsum(tmap)
  
  tLD <- LDshrink::LDshrink(haplomat, mapdat, na.rm = F, m=85, Ne=11490.672741, cutoff = 0.001)
  sLD <- LDshrink::sparse_LDshrink(data = nhaplomat, indices=0:(p-1),total_size=p,
                         mapd = mapdat, m=85, Ne=11490.672741, cutoff = 0.001, useLDshrink=T, progress=F)
  dsLD <- as.matrix(sLD)
  
  summary(c(tLD-dsLD))
  expect_equivalent(tLD, dsLD, 1e-4)
})

test_that("Sparse implementations are equivalent even when chunked", {
  n <- 501
  p <- 1100
  haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
  tmap <- runif(p)
  mapdat <- cumsum(tmap)
  rs <- as.character(1:p)
  shap <- scale(haplomat, center=T, scale=F)
  sLD <- ld2df(data = haplomat,
               mapd = mapdat,
                rsid = as.character(1:p),
               m=85, Ne=11490.672741, cutoff = 0.001, r2cutoff = 0, useLDshrink=T) %>%
    dplyr::arrange(rowsnp, colsnp, r)
  hap_a <- haplomat[, 1:500]
  hap_b <- haplomat[, -(1:500)]
  map_a <- mapdat[1:500]
  map_b <- mapdat[-(1:500)]
  rs_a <- rs[1:500]
  rs_b <- rs[-(1:500)]
  
  df_ab <- ld2df_p(data_a = hap_a,
                   data_b = hap_b,
                   mapd_a = map_a,
                   mapd_b = map_b,
                   rsid_a = rs_a,
                   rsid_b = rs_b,
                   m=85,
                   Ne=11490.672741,
                   cutoff = 0.001, r2cutoff = 0,
                   useLDshrink = T)
  df_aa <- ld2df(data = hap_a,
                 mapd = map_a,
                 rsid = rs_a,
                 m=85,
                 Ne=11490.672741,
                 cutoff = 0.001, r2cutoff = 0,
                 useLDshrink=T)
  df_bb <- ld2df(data = hap_b,
                 mapd = map_b,
                 rsid = rs_b,
                 m=85,
                 Ne=11490.672741,
                 cutoff = 0.001,r2cutoff = 0,
                 useLDshrink=T)
  df_p <- dplyr::bind_rows(df_aa,df_ab,df_bb) %>% dplyr::arrange(rowsnp,colsnp,r)
  expect_equal(df_p,sLD)
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


# test_that("LD2df correctly subsets", {
#   n <- 500
#   p <- 1100
#   p_tot <- (p^2-p)/2+p
#   haplomat <- matrix(sample(0:1,n*p,replace = T),n,p)
#   tmap <- runif(p)
#   
#   mapdat <- cumsum(tmap)
# 
#   tLDdf <- ld2df(LDshrink(haplomat,mapdat,na.rm = F,m = 85,Ne=11490.672741,cutoff=0.001),rsid = as.character(1:p),r2cutoff = 0.001)
#   nhaplomat <- scale(haplomat,center=T,scale=F)
#   nLDdf <- ld2df_p(scaled_data = nhaplomat,mapd = mapdat,rsid =as.character(1:p),m=85,Ne=11490.672741, cutoff = 0.001,r2cutoff = 0.001)
#   
#   nR <- matrix(NA,p,p)
#   nR[cbind(as.integer(tLDdf$rowsnp),as.integer(tLDdf$colsnp))] <- tLDdf$r2
#   cR <- matrix(NA,p,p)
#   cR[cbind(as.integer(nLDdf$rowsnp),as.integer(nLDdf$colsnp))] <- nLDdf$r2
#   testthat::expect_equal(nrow(nLDdf),nrow(tLDdf))
#   testthat::expect_equal(cR[upper.tri(cR)],nR[upper.tri(nR)],1e-4)
# })



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











