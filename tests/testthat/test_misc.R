context("Testing miscellaneous functions")

test_that("dcor", {
  expect_equal(dcor(t(matrix(rep(1:10, 2), ncol = 2)))[1], 0)
})

test_that("getScoreMatrix", {
  m <- getScoreMatrix("BLOSUM62")
  expect_true(nrow(m) == 25L)
  expect_true(ncol(m) == 25L)
  expect_true(m["H", "H"] == 8L)
  expect_true(m["A", "W"] == -3L)
})