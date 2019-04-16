context("checking alleles")


test_that("don't flip alleles for the same dataframe",{
  p <- 100
  sample_df <- tibble::tibble(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(size=p,replace=T) %>% tidyr::unite(allele,sep=",")
  flip_res <- flip_alleles(sample_df$allele,sample_df$allele)
  expect_equal(flip_res,rep(1,nrow(sample_df)))
})
test_that("flip all alleles for the same dataframe (but flipped)",{
  p <- 100
  sample_df <- tibble::tibble(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(size=p,replace=T) 
  snp_df <- sample_df %>%tidyr::unite(allele,sep=",")
  ld_df <- dplyr::select(sample_df,alt,ref) %>%tidyr::unite(allele,sep=",")
  flip_res <- flip_alleles(snp_df$allele,ld_df$allele)
  
  expect_equal(flip_res,rep(-1,nrow(snp_df)))
})

test_that("report 0 when alleles don't match", {
  p <- 100
  sample_df <- tibble::tibble(ref=sample(c("A", "C", "G", "T"), 100, replace=T),
                                  alt=sample(c("A", "C", "G", "T"), 100, replace=T)) %>%
    dplyr::filter(ref!=alt) %>%
    dplyr::sample_n(p, replace=T)
  snp_df <- sample_df %>%
    tidyr::unite(allele, sep=",")
  sample_df$ref[4] <- "F"
  sample_df$alt[6] <- "I"
  ld_df <- dplyr::select(sample_df, alt, ref) %>%
    tidyr::unite(allele, sep=",")
  ret <- flip_alleles(snp_df$allele, ld_df$allele)
  expect_equal(ret[c(4, 6)],c(0,0))
  expect_equal(ret[-c(4, 6)],rep(-1,nrow(snp_df)-2))
})
