context("Testing motif coverage")

f <- system.file("files/meme_ras/meme.xml", package = "motiftools")
x <- readMEME(f)
cov <- getMotifCoverage(x)

test_that("getMotifCoverage works as expected", {
  expect_is(cov, "numeric")
  expect_equal(length(cov), 4)
})

totcov <- getMotifTotalCoverage(x)

test_that("getMotifTotalCoverage works as expected", {
  expect_is(totcov, "numeric")
  expect_equal(length(totcov), 1)
})