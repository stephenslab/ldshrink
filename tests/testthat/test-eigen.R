context("eigen")

test_that("Splitting of LD matrices works when chunksize==p", {
  data("shrink_R")
  p <- ncol(shrink_R)
  chunksize <- p
  big_chunks <- make_chunks(shrink_R,chunksize)
  expect_equal(big_chunks$index[[1]],1:p)
  expect_equal(big_chunks$matrices[[1]],shrink_R)
})


test_that("Splitting of LD matrices works when chunksize==p/2", {
  data("shrink_R")
  p <- ncol(shrink_R)
  chunksize <- p/2
  med_chunks <- make_chunks(shrink_R,chunksize)
  expect_equal(med_chunks$index[[1]],1:chunksize)
  nvec <- (chunksize+1):p
  expect_equal(med_chunks$index[[2]],nvec)
  expect_equal(med_chunks$matrices[[2]],shrink_R[nvec,nvec])
})
