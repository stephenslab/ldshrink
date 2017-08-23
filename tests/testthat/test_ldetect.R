context("LDetect")



test_that("ldetect chunking works like LDshrink chunking",{
  data("map_dat")

  split_ldd <- ldetect_partition(map_dat)
  split_ld <- chunk_LD(map_dat, progress = T)
  expect_equal(split_ld$startpos,split_ldd$startpos)
  expect_equal(split_ld$endpos,split_ldd$endpos)
})
