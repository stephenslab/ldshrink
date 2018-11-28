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




# 
# 
# test_that("can map variants to LD regions",{
# 
#   data("break_df")
#   bad_reg <- dplyr::group_by(break_df,chr) %>% dplyr::slice(c(1)) %>% dplyr::mutate(stop=start,start=0) %>% dplyr::ungroup()
#   bad_reg2 <- dplyr::group_by(break_df,chr) %>% dplyr::slice(n()) %>% dplyr::mutate(start=stop,stop=stop+100) %>% dplyr::ungroup()
#   bad_reg <- dplyr::bind_rows(bad_reg,bad_reg2)
#   make_snps <- function(chr,start,stop,...){
#     n_snps <- sample(0:min(c(10,stop-start)),1)
#     return(tibble::data_frame(chr=rep(chr,n_snps),pos=sort(sample((start+1):(stop-1),n_snps,replace=F))))
#   }
#   bad_snps <- purrr::pmap_dfr(bad_reg,make_snps) %>% dplyr::mutate(isbad=T)
#   good_snps <- purrr::pmap_dfr(break_df,make_snps) %>% dplyr::mutate(isbad=F)
#   all_snps <- dplyr::bind_rows(good_snps,bad_snps) %>% dplyr::arrange(chr,pos) %>% dplyr::distinct(chr,pos,.keep_all = T)
# 
# 
#   all_snps <- assign_snp_block(snp_df = all_snps,break_df = break_df,assign_all=F)
#   # all_snps <- dplyr::mutate(all_snps,region_id=reg_id)
# 
# 
#   check_reg <- dplyr::inner_join(all_snps,break_df)
#   # dplyr::filter(check_reg,!dplyr::between(pos,start,stop))
#   good_res <- purrr::pmap_lgl(check_reg,function(pos,start,stop,...){
#     dplyr::between(pos,start,stop)
#   })
#   expect_equal(all(good_res),T)
# })
# 


test_that("can chunk snp_df by number of chunks or chunksize", {
  
  # p <- 1000
    gen_chrom <- function(ch, p){

    tibble::data_frame(pos=sort(sample(1:(p*100), p, replace = F))) %>% dplyr::mutate(chr=ch)
  }
  mp <- sample(10^(seq(3, 5, length.out = 50)), 22, replace=F)
  snp_df <- purrr::map2_dfr(1:22, mp, gen_chrom)
  chunk_df <- chunk_genome(snp_df, n_chunks = 10)
  all_chunk_df <- dplyr::group_by(chunk_df, chr) %>% dplyr::summarise(n_chunks=dplyr::n_distinct(region_id))
  expect_equal(unique(all_chunk_df$n_chunks), 10)
  chunk_df <- chunk_genome(snp_df, chunk_size = 100, min_size = 10)
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
  testmap <- interpolate_genetic_map(mymap,pos,not_p)
  expect_equal(not_data,testmap)
})



