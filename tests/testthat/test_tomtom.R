context("Testing readTOMTOM")

f <- system.file("files/tomtom.xml", package = "motiftools")
x <- readTOMTOM(f)

test_that("slots", {
  expect_identical(nquery(x), 9L)
  expect_identical(ntarget(x), 9L)
})

test_that("plotMotifMatches works", {
  p <- plotMotifMatches(x)
  expect_is(p, "gg")
  p <- plotMotifMatches(x, color = "black")
  expect_is(p, "gg")
})