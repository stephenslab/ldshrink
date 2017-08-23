context("Diagonal Approximation")

test_that("cutting out diagonals works as expected",{
  
  diagmat <- matrix(as.numeric(1:81),9,9)
  darr <- approx_diag(diagmat,3)
  expect_equal(diagmat[1:3,1:3],darr[,,1])
  expect_equal(diagmat[(4:6),(4:6)],darr[,,2])
  expect_equal(diagmat[(7:9),(7:9)],darr[,,3])

})