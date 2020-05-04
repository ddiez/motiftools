context("Testing readFIMO")

f <- system.file("files/fimo.xml", package = "motiftools")
x <- readFIMO(f)

test_that("slots", {
  expect_identical(nmotif(x), 9L)
  expect_identical(nseq(x), 4L)
  expect_identical(nmotif(x[, 1]), 1L)
  expect_identical(nseq(x[1, ]), 1L)
  expect_equal(sequenceNames(x), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_equal(motifNames(x), as.character(1:9))
})