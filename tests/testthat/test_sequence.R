library(motiftools)
library(Biostrings)

f <- system.file("files/ras.faln", package = "motiftools")
aln <- readAAMultipleAlignment(f)
cons <- conservationMatrix(aln)

test_that("conservation", {
  expect_type(cons, "integer")
  expect_true(cons[1,1] == 1L)
  expect_true(cons[1,14] == 5L)
})