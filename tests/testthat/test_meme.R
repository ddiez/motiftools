library(motiftools)

f <- system.file("files/meme_ras/meme.xml", package = "motiftools")
x <- readMEME(f)

test_that("slots", {
  expect_identical(motiftools:::nmotif(x), 9L)
  expect_identical(motiftools:::nseq(x), 4L)
  expect_equal(motiftools:::sequenceNames(x), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_equal(motiftools:::motifNames(x), as.character(1:9))
})

