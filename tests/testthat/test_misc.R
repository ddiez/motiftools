context("Testing miscellaneous functions")

test_that("dcor", {
  x <- dcor(t(matrix(rep(c(1:10, 10:1), 2), ncol = 4)))
  expect_is(x, "dist")
  expect_equal(x[1], 2)
  expect_equal(x[2], 0)
})

test_that("getScoreMatrix", {
  m <- getScoreMatrix("BLOSUM62")
  expect_true(nrow(m) == 25L)
  expect_true(ncol(m) == 25L)
  expect_true(m["H", "H"] == 8L)
  expect_true(m["A", "W"] == -3L)
})