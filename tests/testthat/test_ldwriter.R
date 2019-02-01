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





test_that("dense to sparse matrix",{
    library(ldshrink)
    library(testthat)
    Rr <- cor(matrix(rnorm(100*10),100,10))
    testR <- ldshrink:::round_trip_skyline_t(Rr,"dsCMatrix")
    sR <- as(Rr,"dsCMatrix")
    expect_equal(sR,testR)

})
