context("ld_writer")

test_that("round trip dense matrix",{
  pa <- 10
  pb <- 10
  n <- 100
  tx <- matrix(rnorm(n*pa),n,pa)
  Rr <- cor(tx)
  testR <- ldshrink:::round_trip_skyline_t(Rr,"matrix")
  expect_equal(testR,Rr,check.attributes=FALSE)
})



test_that("round trip dataframe ",{
  library(ldshrink)
  library(testthat)
  Rr <- cor(matrix(rnorm(100*10),100,10))
  testR <- ldshrink:::round_trip_skyline_t(Rr,"data.frame")
  tdf <- tidyr::unnest(testR)
  expect_equal(Rr[cbind(tdf$rowsnp,tdf$colsnp)],tdf$r)
  newdf <- dplyr::filter(tibble::data_frame(rowsnp=c(row(Rr)),colsnp=c(col(Rr)),r=c(Rr)),rowsnp<=colsnp)
  expect_equal(newdf,tdf)
})


test_that("dense to sparse matrix",{
    library(ldshrink)
    library(testthat)
    Rr <- cor(matrix(rnorm(100*10),100,10))
    testR <- ldshrink:::round_trip_skyline_t(Rr,"dsCMatrix")
    sR <- as(Rr,"dsCMatrix")
    expect_equal(sR,testR)

})
