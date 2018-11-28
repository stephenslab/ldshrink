context("checking alleles")


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
  ld_df <- dplyr::select(sample_df,alt,ref) %>%tidyr::unite(allele,sep=",")
  expect_true(all(flip_allele_exp(snp_df$allele,ld_df$allele)))
})

test_that("report NA when alleles don't match", {
  p <- 100
  sample_df <- tibble::data_frame(ref=sample(c("A", "C", "G", "T"), 100, replace=T),
                                  alt=sample(c("A", "C", "G", "T"), 100, replace=T)) %>%
    dplyr::filter(ref!=alt) %>%
    dplyr::sample_n(p, replace=T)
  snp_df <- sample_df %>%
    tidyr::unite(allele, sep=",")
  sample_df$ref[4] <- "F"
  sample_df$alt[6] <- "I"
  ld_df <- dplyr::select(sample_df, alt, ref) %>%
    tidyr::unite(allele, sep=",")
  ret <- flip_allele_exp(snp_df$allele, ld_df$allele)
  expect_true(all(ret[-c(4, 6)]))
  expect_true(all(is.na(ret[c(4, 6)])))
})
