library(motifTools)

test_that("alignment", {
  library(Biostrings)
  data(BLOSUM62)
  expect_equal(as.character(align("ALA", "ALA", score.matrix = BLOSUM62, type = "global")$alignment[1, 1]), "A")
})