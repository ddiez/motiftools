library(motiftools)
library(Biostrings)

f <- system.file("files/ras.faln", package = "motiftools")
aln <- readAAMultipleAlignment(f)

test_that("conservation", {
  expect_equal(conservationMatrix(aln)[1,1], 1L)
  expect_equal(conservationMatrix(aln)[1,14], 5L)
})