context("utils")

library(dplyr)
data(break_df)



test_that("don't flip alleles for the same dataframe",{
  p <- 100
  sample_df <- tibble::data_frame(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(p,replace=T) %>% tidyr::unite(allele,sep=",")
  
  expect_true(all(!flip_allele_exp(sample_df$allele,sample_df$allele)))
})
test_that("flip all alleles for the same dataframe (but flipped)",{
  p <- 100
  sample_df <- tibble::data_frame(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(p,replace=T) 
  snp_df <- sample_df %>%tidyr::unite(allele,sep=",")
  ld_df <- select(sample_df,alt,ref) %>%tidyr::unite(allele,sep=",")
  expect_true(all(flip_allele_exp(snp_df$allele,ld_df$allele)))
})

test_that("report NA when alleles don't match",{
  p <- 100
  sample_df <- tibble::data_frame(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(p,replace=T) 
  snp_df <- sample_df %>%tidyr::unite(allele,sep=",")
  sample_df$ref[4] <- "F"
  sample_df$alt[6] <- "I"
  ld_df <- dplyr::select(sample_df,alt,ref) %>%tidyr::unite(allele,sep=",")
  ret <- flip_allele_exp(snp_df$allele,ld_df$allele)
  expect_true(all(ret[-c(4,6)]))
  expect_true(all(is.na(ret[c(4,6)])))
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
  
  
  all_snps <- assign_snp_block(snp_df = all_snps,break_df = break_df,assign_all=F)
  # all_snps <- dplyr::mutate(all_snps,region_id=reg_id)

  
  check_reg <- dplyr::inner_join(all_snps,break_df)
  # dplyr::filter(check_reg,!dplyr::between(pos,start,stop))
  good_res <- purrr::pmap_lgl(check_reg,function(pos,start,stop,...){
    dplyr::between(pos,start,stop)
  })
  expect_equal(all(good_res),T)
})





  
  