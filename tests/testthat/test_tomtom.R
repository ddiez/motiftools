library(motiftools)

f <- system.file("files/tomtom_ras/tomtom.xml", package = "motiftools")
x <- readTOMTOM(f)

test_that("slots", {
  expect_identical(nquery(x), 9L)
  expect_identical(ntarget(x), 9L)
})
