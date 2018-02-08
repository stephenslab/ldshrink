context("LD")

test_that("theta is computed correctly",{
  m <- 100
  nmsum = sum(1 / (1:(2*m-1)))
  theta = (1/nmsum) / (2*m + 1/nmsum)
  expect_equal(calc_theta(m),theta)
})


# test_that("LD shrinkage estimators work as expected on simulated data",{
#   m <- 100
#   Ne <- 10000
#   n <- 100
#   p <- 500
#   cutoff <- 1e-3
#   tmap <- cumsum(runif(p)/10)
#   
#   # library(RcppEigenH5) 
#   #  ikgf <- "/media/nwknoblauch/Data/GTEx/GTEx_rssr/Genome_SNP/SNP_Whole_Blood_1kg_6250_t.h5"
#   #  cuda_f <- "/home/nwknoblauch/Dropbox/test_cuda/test_cov.h5"
#   #  cuda_ijk <- read_df_h5(cuda_f,"R",subcols=c("k","j","i"))
#   #  cuda_sighat <- read_2d_mat_h5(cuda_f,"R","LDshrink")
#   #  cuda_cor <- cov2cor(cuda_sighat)
#   #  matA <- read_2d_index_h5(ikgf,"SNPdata","genotype",1:6000)
#   #  mapA <- read_dvec(ikgf,"SNPinfo","map")[1:6000]
#   # m=85
#   # Ne=1490.672741
#   # cutoff=1e-3
#   # nLD <- calcLD(matA,mapA,m,Ne,cutoff)
#   
#   
#   #   
#   
#   Hpanel <- matrix(sample(c(0,1),n*2*p,replace=T),n*2,p)
#   # mfile <- system.file("m_files/run_install.m",package="rssr")
#   mdir <- system.file("m_files",package="RSSReQTL")
#   
#   #change to the directory with the .m files in Octave
#   library(RcppOctave)
#   .CallOctave('cd',mdir)
#   msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
#   
#   Rmsig <- cov2cor(msig)
#   # Rmsig[lower.tri(Rmsig)] <- 0
#   Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
#   
#   which(abs((Rsig-Rmsig))==max(abs(Rsig-Rmsig)),arr.ind = T)
#   
#   # Rsig[lower.tri(Rsig)] <- 0
#   expect_equal(Rsig,Rmsig)
# })

test_that("LD shrinkage estimators results are approximately equal betweeen haplotype and genotype data",{
  m <- 100
  Ne <- 10000
  n <- 100
  p <- 500
  cutoff <- 1e-3
  tmap <- cumsum(runif(p)/10)
  
  # library(RcppEigenH5) 
  #  ikgf <- "/media/nwknoblauch/Data/GTEx/GTEx_rssr/Genome_SNP/SNP_Whole_Blood_1kg_6250_t.h5"
  #  cuda_f <- "/home/nwknoblauch/Dropbox/test_cuda/test_cov.h5"
  #  cuda_ijk <- read_df_h5(cuda_f,"R",subcols=c("k","j","i"))
  #  cuda_sighat <- read_2d_mat_h5(cuda_f,"R","LDshrink")
  #  cuda_cor <- cov2cor(cuda_sighat)
  #  matA <- read_2d_index_h5(ikgf,"SNPdata","genotype",1:6000)
  #  mapA <- read_dvec(ikgf,"SNPinfo","map")[1:6000]
  # m=85
  # Ne=1490.672741
  # cutoff=1e-3
  # nLD <- calcLD(matA,mapA,m,Ne,cutoff)
  
  
  #   
  
  Gpanel <- matrix(sample(c(0,2),n*2*p,replace=T),n*2,p)
  Hpanel <- Gpanel/2
  RHsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  RGsig <- calcLD(hmata = Gpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  
  # Rsig[lower.tri(Rsig)] <- 0
  expect_equal(RGsig,RHsig,tolerance=1e-2)
})



# 
# test_that("LD shrinkage estimators work the same real data ",{
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
#   mdir <- system.file("m_files",package="RSSReQTL")  
#   #change to the directory with the .m files in Octave
#   library(RcppOctave)
#   .CallOctave('cd',mdir)
#   msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
#   Rmsig <- cov2cor(msig)
#   Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
#   evdRsig <- eigen(Rsig)
#   evdmsig <- eigen(Rmsig)
#   min(evdmsig$values)
#   min(evdRsig$values)
#   which(abs(Rsig-Rmsig)==max(abs(Rsig-Rmsig)),arr.ind = T)
#     expect_equal(Rsig,Rmsig)
# })




# 
# test_that("LD shrinkage estimators work the same real (larger) data ",{
#   
#   m=85
#   Ne=1490.672741
#   cutoff=1e-3
#   data("bighaplo")
#   data("bigmap")
#   Hpanel <- bighaplo
#   tmap <- bigmap
#   
#   
#   mdir <- system.file("m_files",package="RSSReQTL")  
#   #change to the directory with the .m files in Octave
#   library(RcppOctave)
#   .CallOctave('cd',mdir)
#   msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
#   Rmsig <- cov2cor(msig)
#   
#   Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
#   
#   expect_equal(Rsig,Rmsig)
# })

# 
# test_that("LD shrinkage estimators work the same real (larger) data and a big cutoff ",{
#   
#   m=85
#   Ne=1490.672741
# 
#   data("bighaplo")
#   data("bigmap")
#   Hpanel <- bighaplo
#   tmap <- bigmap
#   
#   
#   mdir <- system.file("m_files",package="RSSReQTL")  
#   #change to the directory with the .m files in Octave
#   library(RcppOctave)
#   .CallOctave('cd',mdir)
#   cutoff=.5
#   msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
#   Rmsig <- cov2cor(msig)
#   
#   Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
#   
#   expect_equal(Rsig,Rmsig)
#   
#   
# })





test_that("LD shrinkage estimators give similar results for sparse and dense data",{
  
  m=85
  Ne=11490.672741
  cutoff=1e-3
  data("mapdat")
  data("haplomat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  Rsig_h_d <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig_h_s <- as.matrix(sp_calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff))
  expect_equivalent(Rsig_h_d,Rsig_h_s)
  
})

# test_that("Distance can be computed using only 'linear algebra'",{
#   
#   data("mapdat")
#   tmap <- mapdat  
#   tmap <- t(t(tmap))
#   Rdist <- abs(outer(tmap,tmap,"-"))
#   onemat <- matrix(-1,length(tmap),length(tmap))
#   ndist <- t(tmap)%*%onemat%*%tmap
#   system.time({
#   tlf <- "~/Desktop/scratch/polyg_scratch/h5/sub_seq_hapmap_haplo.h5"
#   hmata <- t(SeqSupport::read_2d_mat_h5(tlf,"/","dosage"))
#   map <- SeqSupport::read_vec(tlf,"SNPinfo/map")
#   tLD <- calcLD(hmata = hmata,mapa = map)
#   SeqSupport::write_mat_h5("/media/nwknoblauch/Data/testR.h5",groupname = "LD",dataname = "R",data = tLD,deflate_level = 4)
#   })
# })




