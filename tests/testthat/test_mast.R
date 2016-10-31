library(motiftools)

f <- system.file("files/mast_ras/mast.xml", package = "motiftools")
x <- readMAST(f)

test_that("slots", {
  expect_identical(nmotif(x), 9L)
  expect_identical(nseq(x), 4L)
  expect_identical(nmotif(x[, 1]), 1L)
  expect_identical(nseq(x[1, ]), 1L)
  expect_equal(sequenceNames(x), c("RASK_MOUSE", "RASN_MOUSE", "RASH_MOUSE", "RASM_HUMAN")) # MAST gives a different order (ordered by evalue).
  expect_equal(motifNames(x), as.character(1:9))
})

