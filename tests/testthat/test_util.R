context("utils")

library(dplyr)
data(break_df)

test_that("I can calculate correlation with MKL",{
  
  n <- 100
  p <- 1000
  tmat <- matrix(rnorm(n*p),n,p)
  retC <- cov_mkl(tmat)
  rC <- cov(tmat)
  expect_equal(rC,retC)
  
})


test_that("Two equivalent ways to calculate covariance",{
  
  n <- 100
  p <- 1000
  tmat <- matrix(rnorm(n*p),n,p)
  retC <- calc_cov(tmat)
  retC2 <- calc_cov_s(tmat)
  expect_equal(retC,retC2)
  
  res <- microbenchmark(r=calc_cov(tmat),
                 s=calc_cov_s(tmat))
  
})

test_that("can break up LD regions using ranges",{
  
  ldr <- 1:10
  ldv <- integer()
  for(i in ldr){
    ldv <- c(ldv,rep(i,sample(0:10,1)))
  }
  res <- ld_chunks(ldv)
  r_res <- split(ldv,ldv) %>% map_int(length)
  
  expect_equivalent(res,r_res)
  
})

test_that("can map variants to LD regions",{

  data("break_df")
  bad_reg <- dplyr::group_by(break_df,chr) %>% dplyr::slice(c(1)) %>% dplyr::mutate(stop=start,start=0) %>% dplyr::ungroup()
  bad_reg2 <- dplyr::group_by(break_df,chr) %>% dplyr::slice(n()) %>% dplyr::mutate(start=stop,stop=stop+100) %>% dplyr::ungroup()
  bad_reg <- dplyr::bind_rows(bad_reg,bad_reg2)
  make_snps <- function(chr,start,stop,...){
    n_snps <- sample(0:min(c(10,stop-start)),1)
    return(tibble::data_frame(chr=rep(chr,n_snps),pos=sort(sample((start+1):(stop-1),n_snps,replace=F))))
  }
  bad_snps <- purrr::pmap_dfr(bad_reg,make_snps) %>% dplyr::mutate(isbad=T)
  good_snps <- purrr::pmap_dfr(break_df,make_snps) %>% dplyr::mutate(isbad=F)
  all_snps <- dplyr::bind_rows(good_snps,bad_snps) %>% dplyr::arrange(chr,pos) %>% dplyr::distinct(chr,pos,.keep_all = T)
  
  all_reg <- dplyr::bind_rows()
  reg_id <- set_ld_region(ld_regions = break_df,snp_info = all_snps)
  all_snps <- dplyr::mutate(all_snps,region_id=reg_id)

  
  check_reg <- dplyr::inner_join(all_snps,break_df)
  dplyr::filter(check_reg,!dplyr::between(pos,start,stop))
  good_res <- purrr::pmap_lgl(check_reg,function(pos,start,stop,...){
    dplyr::between(pos,start,stop)
  })
  expect_equal(all(good_res),T)
})





  
  