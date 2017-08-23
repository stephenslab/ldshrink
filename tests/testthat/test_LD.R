context("Test LD calculations")


test_that("Covariance with two input matrices works as expected",{
  
  inp_a <- matrix(runif(90*10),90,10)
  inp_b <- matrix(runif(90*9),90,9)  
  
  R_cov <- cov(inp_a,inp_b)
  
  ld_cov <- calc_cov_p(inp_a,inp_b)
  expect_equal(ld_cov,R_cov)
  
})

test_that("Trivial distributed case works",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  tRsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  theta <- calc_theta(m = m)
  ldp <- c(m,Ne,cutoff,theta)
  p <- ncol(Hpanel)
  chunkp <- as.integer(c(0,0,p))
  dRsig <- calcLD_par(hmat = Hpanel,map = tmap,ldparams = ldp,id = chunkp)

  expect_equal(tRsig,dRsig)
})


test_that("off-diagonal LD blocks can be calculated",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  tRsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  theta <- calc_theta(m = m)
  ldp <- c(m,Ne,cutoff,theta)
  p <- ncol(Hpanel)
  csize <- 501
  chunkp <- as.integer(c(0L,1L,csize))
  dRsig <- calcLD_par(hmat = Hpanel,map = tmap,ldparams = ldp,id = chunkp)
  expect_true(all(abs(dRsig)<=1))

  # Rsig[lower.tri(Rsig)] <- 0
  sub_tRsig <- tRsig[1:501,502:1000]
  expect_equal(sub_tRsig,dRsig)
  
})


test_that("diagonal LD blocks can be calculated",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  tRsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  theta <- calc_theta(m = m)
  ldp <- c(m,Ne,cutoff,theta)
  p <- ncol(Hpanel)
  csize <- as.integer(p/2)
  chunkp <- as.integer(c(1L,1L,csize))
  dRsig <- calcLD_par(hmat = Hpanel,map = tmap,ldparams = ldp,id = chunkp)
  # Rsig[lower.tri(Rsig)] <- 0
  
  expect_equal(tRsig[501:1000,501:1000],dRsig)
  
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


test_that("MPI works as expected",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  tRsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rr <- tRsig[1:10,1:10]
  mpir <- matrix(scan("inst/mpi_example/LDmat.txt",what=numeric()),10,10,byrow = T)
  expect_equal(Rr,mpir)
  
})








