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



test_that("ld scores are calculated correctly",{
  data("reference_genotype")  
  data("reference_map")
  data("reference_ldscores")
  RLD <- estimate_LD(reference_panel = reference_genotype,map = reference_map,method = "sample")
  L2 <- estimate_LDscores(RLD,nrow(reference_genotype))
  expect_equal(L2,reference_ldscores,check.attributes=F,tolerance=1e-3)
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

test_that("genetic map interpolation works when query is a subset ",{
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

test_that("genetic map interpolation works when query ==target",{
  p <- 10L
  b <- runif(1)
  pos <- 1L:p
  not_pos <- pos
  mymap <- pos*b
  testmap <- interpolate_genetic_map(mymap,pos,not_pos,strict=T)
  expect_equal(testmap,mymap)
})


test_that("Check for assigning SNPs to blocks",{
  
  input_f <- system.file("test_data/fourier_ls-all.bed.gz",package = "ldshrink")
  ld_df <- read_tsv(input_f,col_types=cols(
    chr = col_character(),
    start = col_integer(),
    stop = col_integer()
  )) %>%
    group_by(chr) %>%
    mutate(start = if_else(start==min(start),
                           0L,start)) %>%
    mutate(stop = if_else(stop==max(stop),
                          .Machine$integer.max,stop)) %>%
    ungroup() %>%
    mutate(region_id=1:n())
  

  snp_df <- pmap_dfr(ld_df,function(chr,start,stop,region_id,...){
    n_snps <- sample(1:100,1)
    tibble(pos=sort(sample(start:stop,n_snps,replace=T))) %>% 
      mutate(chr=chr,region_id=region_id) 
  })
  n_region_id <- assign_region(break_chr = ld_df$chr,
                               break_start = ld_df$start,
                               break_stop = ld_df$stop,
                               break_id = ld_df$region_id,
                               snp_chr = snp_df$chr,
                               snp_pos = snp_df$pos,assign_all=T)
                               
                               # break_df = ld_df,assign_all = TRUE)
  expect_equal(snp_df$region_id,n_region_id)
})

test_that("Check for allele flipping works",{
  
  nuc <- c("A","C","T","G")
  ref_a <- outer(nuc,nuc,function(x,y){paste0(x,",",y)})
  ref_a <- ref_a[ref_a!=diag(ref_a)]
  
  #don't flip yourself
  expect_equal(flip_alleles(ref_a,ref_a),rep(1L,12))
  rev_ref <- purrr::map_chr(strsplit(ref_a,",",fixed=T),~paste(rev(.x),collapse=",")) 
  expect_equal(flip_alleles(ref_a,rev_ref),rep(-1L,12))
  
  flips <- c("A"="T",
             "T"="A",
             "G"="C",
             "C"="G")
  
  flip_ref <- purrr::map_chr(strsplit(ref_a,",",fixed=T),~paste(flips[.x],collapse=",")) 
  expect_equal(flip_alleles(ref_a,flip_ref),rep(-1L,12))
  flip_rev_ref <- purrr::map_chr(strsplit(ref_a,",",fixed=T),~paste(rev(flips[.x]),collapse=",")) 
  expect_equal(flip_alleles(ref_a,flip_rev_ref),rep(1L,12))
  
})



