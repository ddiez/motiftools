context("Testing miscellaneous functions")

test_that("dcor", {
  expect_equal(dcor(t(matrix(rep(1:10, 2), ncol = 2)))[1], 0)
})