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


test_that("don't flip alleles for the same dataframe",{
  p <- 100
  sample_df <- tibble::tibble(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(p,replace=T) %>% tidyr::unite(allele,sep=",")
  
  expect_true(all(!flip_allele_exp(sample_df$allele,sample_df$allele)))
})
test_that("flip all alleles for the same dataframe (but flipped)",{
  p <- 100
  sample_df <- tibble::tibble(ref=sample(c("A","C","G","T"),100,replace=T),alt=sample(c("A","C","G","T"),100,replace=T)) %>% 
    dplyr::filter(ref!=alt) %>% dplyr::sample_n(p,replace=T) 
  snp_df <- sample_df %>%tidyr::unite(allele,sep=",")
  ld_df <- dplyr::select(sample_df,alt,ref) %>%tidyr::unite(allele,sep=",")
  expect_true(all(flip_allele_exp(snp_df$allele,ld_df$allele)))
})

test_that("report NA when alleles don't match", {
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
  ret <- flip_allele_exp(snp_df$allele, ld_df$allele)
  expect_true(all(ret[-c(4, 6)]))
  expect_true(all(is.na(ret[c(4, 6)])))
})






