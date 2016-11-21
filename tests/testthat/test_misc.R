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

test_that("eud works as expected", {
  x <- eud(0, 0)
  expect_is(x, "numeric")
  expect_equal(x, 0)
  expect_equal(eud(0, 1), 1)
  expect_equal(c(0, 1), c(1, 0), sqrt(2))
})