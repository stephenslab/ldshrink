context("utils")


test_that("linear indexing of symmetric matrices",{
  
  p <- 9
  ptot <- (p^2-p)/2+p
  mymat <- matrix(0,p,p)
  mymat[lower.tri(mymat,diag=T)] <- 1:ptot
  
  indfun <- function(index,p){
    k <- index-1
    i = floor( ( 2*p+1 - sqrt( (2*p+1)*(2*p+1) - 8*k ) ) / 2 )
    # i = floor( ( 2*p+1 - sqrt( (2*p+1)*(2*p+1) - 8*k ) ) / 2 ) 
    # j = k - p*i + i*(i-1)/2 
    j = k - (2*p-1- i)*i/2
    return(matrix(c(j,i),1,2)+1)
  }
  for(i in 1:ptot){
    ij <- indfun(i,p)
    expect_equal(mymat[ij],i)
  }
  
  
})


test_that("can chunk snp_df by number of chunks or chunksize", {
    gen_chrom <- function(ch, p){
    tibble::data_frame(pos=sort(sample(1:(p*100), p, replace = F))) %>% dplyr::mutate(chr=ch)
  }
  mp <- sample(10^(seq(3, 5, length.out = 50)), 22, replace=F)
  snp_df <- purrr::map2_dfr(1:22, mp, gen_chrom)
  chunk_df <- chunk_genome(snp_df, n_chunks = 10)
  all_chunk_df <- dplyr::group_by(chunk_df, chr) %>% dplyr::summarise(n_chunks=dplyr::n_distinct(region_id))
  expect_equal(unique(all_chunk_df$n_chunks), 10)
  chunk_df <- chunk_genome(snp_df, chunk_size = 100,min_size = 10)
  all_chunk_df <- dplyr::group_by(chunk_df, region_id) %>% dplyr::summarise(region_size=n())
  expect_gt(min(all_chunk_df$region_size), 10)
  expect_lte(max(all_chunk_df$region_size), 110)
})

test_that("genetic map interpolation works",{
  p <- 10000
  b <- runif(1)
  full_p <- 1:(100*p)
  pos <- sort(sample(full_p,p,replace=F))
  not_p <- full_p[!full_p %in% pos]
  mymap <- pos*b
  not_data <- not_p*b
  testmap <- interpolate_genetic_map(mymap,pos,not_p,strict=F)
  expect_equal(not_data,testmap)
})


test_that("ldshrink can work like base R for sample correlation LD scores",{
  # n <- 500
  # p <- 1100
  library(ldshrink)
  data("reference_genotype")
  data("reference_map")
  ldscoref <- system.file("reference_genotype.l2.ldscore.gz",package = "ldshrink")
  ld_df <- readr::read_tsv(ldscoref)
  tR <- cor(reference_genotype)
  rResult <- estimate_LD(reference_panel = reference_genotype,method="sample",map = reference_map,output = "matrix")
  L2 <- estimate_LDscores(rResult,nrow(reference_genotype))
  # R <- cor(reference_genotype[,1:5])
  expect_equal(ld_df$L2,L2,check.attributes=F,tolerance=1e-3)
  
})



